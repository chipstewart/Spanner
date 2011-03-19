// ***************************************************************************
// MosaikAlignment.cpp - a concatenation of MOSAIK classes that allows the user
//                     to read from and write to MOSAIK alignment archives.
// ---------------------------------------------------------------------------
// (c) 2006 - 2009 Michael Strömberg
// Marth Lab, Department of Biology, Boston College
// ---------------------------------------------------------------------------
// Dual licenced under the GNU General Public License 2.0+ license or as
// a commercial license with the Marth Lab.
// ***************************************************************************

#include "MosaikAlignment.h"

namespace Mosaik {

	// ================
	// CAlignmentReader
	// ================

	// constructor
	CAlignmentReader::CAlignmentReader(void)
		: mIsOpen(false)
		, mInStream(NULL)
		, mNumReads(0)
		, mNumBases(0)
		, mCurrentRead(0)
		, mReadsOffset(0)
		, mReferenceGapOffset(0)
		, mIndexOffset(0)
		, mBuffer(NULL)
		, mBufferPtr(NULL)
		, mBufferLen(0)
		, mCompressionBuffer(NULL)
		, mCompressionBufferLen(0)
		, mPartitionSize(0)
		, mPartitionMembers(0)
		, mRefSeqLUT(NULL)
		, mStatus(AS_UNKNOWN)
		, mSeqTech(ST_UNKNOWN)
	{}

	// destructor
	CAlignmentReader::~CAlignmentReader(void) {
		if(mIsOpen)            Close();
		if(mBuffer)            delete [] mBuffer;
		if(mCompressionBuffer) delete [] mCompressionBuffer;

		// delete the reference sequence LUT
		for(unsigned short i = 0; i < mNumRefSeqs; ++i) delete [] mRefSeqLUT[i];
		delete [] mRefSeqLUT;
	}

	// checks to see if this is truly an MOSAIK alignment archive
	bool CAlignmentReader::CheckFile(const string& filename, const bool showError) {

		// read in the first 6 characters
		char signature[7];
		signature[6] = 0;
		bool foundError = false;

		const char* MOSAIK_SIGNATURE = "MSKAA\4";

		// open the MOSAIK alignment archive
		FILE* checkStream = fopen(filename.c_str(), "rb");
		if(!checkStream) {
			if(showError) {
				printf("ERROR: Could not open %s when validating the alignment archive.\n", filename.c_str());
				exit(1);
			}

			foundError = true;
		}

		// retrieve the MOSAIK alignment archive signature
		if(!foundError) {

			// check if we were able to read 6 bytes
			if(fread(signature, 1, 6, checkStream) < 6) {
				if(showError) {
					printf("ERROR: It seems that the input file (%s) is not in the MOSAIK read format.\n", filename.c_str());
					exit(1);
				}

				foundError = true;
			}

			// check if the read signatures match
			if(!foundError && (strncmp(signature, MOSAIK_SIGNATURE, 5) != 0)) {
				if(showError) {
					printf("ERROR: It seems that the input file (%s) is not in the MOSAIK alignment format.\n", filename.c_str());
					exit(1);
				}

				foundError = true;
			}

			// check if the file format is from another version
			if(!foundError && (MOSAIK_SIGNATURE[5] != signature[5])) {
				if(showError) {
					printf("ERROR: It seems that the input file (%s) was created in another version of MosaikAligner. A new alignment archive is required.\n", filename.c_str());
					printf("       file version: %hu, expected version: %hu\n", signature[5], MOSAIK_SIGNATURE[5]);
					exit(1);
				}

				foundError = true;
			}
		}

		// close the file
		if(checkStream) fclose(checkStream);

		// return the appropriate values
		if(foundError) return false;
		return true;
	}

	// closes the alignment archive
	void CAlignmentReader::Close(void) {
		mIsOpen = false;
		fclose(mInStream);
	}

	// returns the a pointer to the header tags map
	//map<unsigned char, Tag>* CAlignmentReader::GetHeaderTags(void) {
	//	return &mHeaderTags;
	//}

	// returns the number of reads in the archive
	uint64_t CAlignmentReader::GetNumBases(void) const {
		if(!mIsOpen) return 0;
		return mNumBases;
	}

	// returns the number of reads in the archive
	uint64_t CAlignmentReader::GetNumReads(void) const {
		if(!mIsOpen) return 0;
		return mNumReads;
	}

	// retrieves the read group given a read group code
	ReadGroup CAlignmentReader::GetReadGroupFromCode(const unsigned int code) {
		map<unsigned int, ReadGroup>::const_iterator rgIter = mReadGroupLUT.find(code);

		if(rgIter == mReadGroupLUT.end()) {
			printf("ERROR: The following read group ID was not found in the lookup table: %u\n", code);
			exit(1);
		}

		return rgIter->second;
	}

	// retrieves the read groups vector
	ReadGroupVector CAlignmentReader::GetReadGroups(void) const {
		return mReadGroups;
	}

	// retrieves the reference sequence data
	RefVector CAlignmentReader::GetReferenceData(void) const {
		return mReferenceSequences;
	}

	// retrieves the file status
	AlignmentStatus CAlignmentReader::GetStatus(void) const {
		return mStatus;
	}

	// jumps to the block containing the specified reference index and position
	void CAlignmentReader::Jump(const unsigned int referenceIndex, const unsigned int referencePosition) {

		// ===============
		// parse the index
		// ===============

		if(mIndexOffset == 0) {
			cout << "ERROR: Cannot jump to the desired compressed block because the index offset was not set." << endl;
			exit(1);
		}

		fseek64(mInStream, mIndexOffset, SEEK_SET);

		// read the number of entries
		unsigned int numIndexEntries = 0;
		fread((char*)&numIndexEntries, SIZEOF_INT, 1, mInStream);

		// load the index
		const unsigned int requestedBytes = numIndexEntries * (SIZEOF_UINT64 + SIZEOF_INT + SIZEOF_SHORT);
		CMemoryUtilities::CheckBufferSize(mBuffer, mBufferLen, requestedBytes);
		fread(mBuffer, requestedBytes, 1, mInStream);

		// find the block containing the specified reference index and position
		unsigned int bufferOffset = 0;

		unsigned short index      = 0;
		unsigned int position     = 0;
		off_type offset           = 0;

		bool foundBlock = false;

		for(unsigned int i = 0; i < numIndexEntries; ++i) {

			// retrieve the reference index
			memcpy((char*)&index, mBuffer + bufferOffset, SIZEOF_SHORT);
			bufferOffset += SIZEOF_SHORT;

			// store the reference position
			memcpy((char*)&position, mBuffer + bufferOffset, SIZEOF_INT);
			bufferOffset += SIZEOF_INT;

			// store the file offset
			memcpy((char*)&offset, mBuffer + bufferOffset, SIZEOF_UINT64);
			bufferOffset += SIZEOF_UINT64;

			// keep going until we find a compression block that is past our desired index and position
			if(index > referenceIndex) foundBlock = true;
			if((index == referenceIndex) && (position >= referencePosition)) foundBlock = true;
			if(foundBlock) break;
		}

		if(!foundBlock) {
			cout << "ERROR: A suitable compression block was not found in the index." << endl;
			exit(1);
		}

		fseek64(mInStream, offset, SEEK_SET);
		mCurrentRead      = 0;
		mPartitionMembers = 0;
		mPartitionSize    = 0;
	}

	// loads the next read from the alignment archive
	bool CAlignmentReader::LoadNextRead(Mosaik::AlignedRead& ar) {

		if(!mIsOpen) {
			cout << "ERROR: An attempt was made to get reads from an alignment archive that hasn't been opened yet." << endl;
			exit(1);
		}

		// check if we have already processed all of the reads
		if(mCurrentRead >= mNumReads) return false;

		// read the partition
		if(mPartitionMembers == mPartitionSize) {
			if(!ReadPartition()) return false;
		}

		// initialize
		unsigned char readStatus = RF_UNKNOWN;

		unsigned int numMate1Alignments = 0;
		unsigned int numMate2Alignments = 0;

		// load the read header
		LoadReadHeader(ar.Name, ar.ReadGroupCode, readStatus, numMate1Alignments, numMate2Alignments);

		// interpret the read status
		const bool haveMate1        = ((readStatus & RF_HAVE_MATE1)              != 0 ? true : false);
		const bool haveMate2        = ((readStatus & RF_HAVE_MATE2)              != 0 ? true : false);
		const bool isResolvedAsPair = ((readStatus & RF_RESOLVED_AS_PAIR)        != 0 ? true : false);
		ar.IsLongRead               = ((readStatus & RF_IS_LONG_READ)            != 0 ? true : false);
		ar.IsPairedEnd              = ((readStatus & RF_IS_PAIRED_IN_SEQUENCING) != 0 ? true : false);

		// =================================
		// deserialize each mate 1 alignment
		// =================================

		ar.Mate1Alignments.resize(numMate1Alignments);
		if(haveMate1) ReadAlignments(ar.Mate1Alignments, ar.IsLongRead, ar.IsPairedEnd, isResolvedAsPair, ar.ReadGroupCode);

		// =================================
		// deserialize each mate 2 alignment
		// =================================

		ar.Mate2Alignments.resize(numMate2Alignments);
		if(haveMate2) ReadAlignments(ar.Mate2Alignments, ar.IsLongRead, ar.IsPairedEnd, isResolvedAsPair, ar.ReadGroupCode);

		// increment the read counter
		++mCurrentRead;
		++mPartitionMembers;

		return true;
	}

	// load the read header from disk
	void CAlignmentReader::LoadReadHeader(string& readName, unsigned int& readGroupCode, unsigned char& readStatus, unsigned int& numMate1Alignments, unsigned int& numMate2Alignments) {

		// get the read name
		const unsigned char readNameLength = (unsigned char)*mBufferPtr;
		++mBufferPtr;

		readName.resize(readNameLength);
		memcpy((void*)readName.data(), mBufferPtr, readNameLength);
		mBufferPtr += readNameLength;

		// get the read group code
		memcpy((char*)&readGroupCode, mBufferPtr, SIZEOF_INT);
		mBufferPtr += SIZEOF_INT;

		// get the read status
		readStatus = (unsigned char)*mBufferPtr;
		++mBufferPtr;

		const bool haveMate1 = ((readStatus & RF_HAVE_MATE1) != 0 ? true : false);
		const bool haveMate2 = ((readStatus & RF_HAVE_MATE2) != 0 ? true : false);

		// get the number of mate 1 alignments
		if(haveMate1) {
			memcpy((char*)&numMate1Alignments, mBufferPtr, SIZEOF_INT);
			mBufferPtr += SIZEOF_INT;
		}

		// get the number of mate 2 alignments
		if(haveMate2) {
			memcpy((char*)&numMate2Alignments, mBufferPtr, SIZEOF_INT);
			mBufferPtr += SIZEOF_INT;
		}
	}

	// opens the alignment archive
	void CAlignmentReader::Open(const string& filename) {

		if(mIsOpen) {
			cout << "ERROR: An attempt was made to open an already open alignment archive." << endl;
			exit(1);
		}

		mInputFilename = filename;

		mInStream = fopen(filename.c_str(), "rb");
		if(!mInStream) {
			cout << "ERROR: Could not open the compressed alignment archive (" << mInputFilename << ") for reading." << endl;
			exit(1);
		}

		mIsOpen = true;

		// ===============
		// read the header
		// ===============

		// MOSAIK_SIGNATURE[6]	   0  -  5
		// STATUS[1]               6  -  6
		// SEQUENCE_TECHNOLOGY[2]  7  -  8
		// ARCHIVE_DATE[8]		   9  - 16
		// NUM_REFERENCE_SEQS[4]   17 - 20
		// NUM_READ_GROUPS[4]      21 - 24
		// NUM_READS[8]            25 - 32
		// NUM_BASES[8]            33 - 40
		// REFERENCES_OFFSET[8]    41 - 48
		// REFERENCE_GAP_OFFSET[8] 49 - 57
		// INDEX_OFFSET[8]         58 - 63
		// NUM_READ_GROUP_TAGS[1]  64 - 64
		// READ_GROUPS[*]

		// skip the MOSAIK signature
		const unsigned char SIGNATURE_LENGTH = 6;
		fseek64(mInStream, SIGNATURE_LENGTH, SEEK_SET);

		// retrieve the alignment file status
		mStatus = (AlignmentStatus)fgetc(mInStream);

		// retrieve the sequencing technology
		fread((char*)&mSeqTech, SIZEOF_SHORT, 1, mInStream);

		// skip the archive date
		fseek64(mInStream, SIZEOF_UINT64, SEEK_CUR);

		// retrieve the number of reference sequences
		fread((char*)&mNumRefSeqs, SIZEOF_INT, 1, mInStream);

		// retrieve the number of read groups
		unsigned int numReadGroups;
		fread((char*)&numReadGroups, SIZEOF_INT, 1, mInStream);

		// retrieve the number of reads
		fread((char*)&mNumReads, SIZEOF_UINT64, 1, mInStream);

		// retrieve the number of bases
		fread((char*)&mNumBases, SIZEOF_UINT64, 1, mInStream);

		// retrieve the references offset
		off_type referencesOffset = 0;
		fread((char*)&referencesOffset, SIZEOF_OFF_TYPE, 1, mInStream);

		// retrieve the reference gaps offset
		fread((char*)&mReferenceGapOffset, SIZEOF_OFF_TYPE, 1, mInStream);

		// retrieve the index offset
		fread((char*)&mIndexOffset, SIZEOF_OFF_TYPE, 1, mInStream);

		// retrieve the number of header tags
		const unsigned char numHeaderTags = (unsigned char)fgetc(mInStream);

		if(numHeaderTags != 0) {
			for(unsigned char j = 0; j < numHeaderTags; j++) {
				Tag tag;
				ReadTag(tag);
				mHeaderTags[tag.ID] = tag;
			}
		}

		// DEBUG
		//cout << "mStatus:             " << (short)mStatus << endl;
		//cout << "mSeqTech:            " << mSeqTech << endl;
		//cout << "mNumRefSeqs:         " << mNumRefSeqs << endl;
		//cout << "numReadGroups:       " << numReadGroups << endl;
		//cout << "mNumReads:           " << mNumReads << endl;
		//cout << "mNumBases:           " << mNumBases << endl;
		//cout << "referencesOffset:    " << referencesOffset << endl;
		//cout << "mReferenceGapOffset: " << mReferenceGapOffset << endl;
		//cout << "mIndexOffset:        " << mIndexOffset << endl;
		//cout << "numHeaderTags:       " << (unsigned short)numHeaderTags << endl << endl;

		// retrieve the read groups
		mReadGroups.resize(numReadGroups);

		vector<ReadGroup>::iterator rgIter;
		for(rgIter = mReadGroups.begin(); rgIter != mReadGroups.end(); ++rgIter) {

			// read the metadata string lengths
			const unsigned char centerNameLen   = (unsigned char)fgetc(mInStream);
			const unsigned char libraryNameLen  = (unsigned char)fgetc(mInStream);
			const unsigned char platformUnitLen = (unsigned char)fgetc(mInStream);
			const unsigned char readGroupIDLen  = (unsigned char)fgetc(mInStream);
			const unsigned char sampleNameLen   = (unsigned char)fgetc(mInStream);

			unsigned short descriptionLen = 0;
			fread((char*)&descriptionLen, SIZEOF_SHORT, 1, mInStream);
			fread((char*)&rgIter->SequencingTechnology, SIZEOF_SHORT, 1, mInStream);
			fread((char*)&rgIter->MedianFragmentLength, SIZEOF_INT, 1, mInStream);

			rgIter->CenterName.resize(centerNameLen);
			rgIter->LibraryName.resize(libraryNameLen);
			rgIter->PlatformUnit.resize(platformUnitLen);
			rgIter->ReadGroupID.resize(readGroupIDLen);
			rgIter->SampleName.resize(sampleNameLen);
			rgIter->Description.resize(descriptionLen);

			// read the metadata strings
			fread((void*)rgIter->CenterName.data(),   centerNameLen,   1, mInStream);
			fread((void*)rgIter->Description.data(),  descriptionLen,  1, mInStream);
			fread((void*)rgIter->LibraryName.data(),  libraryNameLen,  1, mInStream);
			fread((void*)rgIter->PlatformUnit.data(), platformUnitLen, 1, mInStream);
			fread((void*)rgIter->ReadGroupID.data(),  readGroupIDLen,  1, mInStream);
			fread((void*)rgIter->SampleName.data(),   sampleNameLen,   1, mInStream);

			// set the read group code
			rgIter->ReadGroupCode = ReadGroup::GetCode(*rgIter);

			// add the read group to our LUT
			mReadGroupLUT[rgIter->ReadGroupCode] = *rgIter;

			// retrieve the number of read group tags
			const unsigned char numReadGroupTags = (unsigned char)fgetc(mInStream);

			if(numReadGroupTags != 0) {
				printf("ERROR: Found %u read group tags, but support for read group tags has not been implemented yet.\n", numReadGroupTags);
				exit(1);
			}

			//// DEBUG
			//cout << "center name:            " << rgIter->CenterName << endl;
			//cout << "description:            " << rgIter->Description << endl;
			//cout << "library name:           " << rgIter->LibraryName << endl;
			//cout << "platform unit:          " << rgIter->PlatformUnit << endl;
			//cout << "read group ID:          " << rgIter->ReadGroupID << endl;
			//cout << "sample name:            " << rgIter->SampleName << endl;
			//cout << "sequencing technology:  " << rgIter->SequencingTechnology << endl;
			//cout << "median fragment length: " << rgIter->MedianFragmentLength << endl << endl;
		}

		// store the reads offset
		mReadsOffset = ftell64(mInStream);

		// ============================
		// read the reference sequences
		// ============================

		// jump to the reference sequence section
		fseek64(mInStream, referencesOffset, SEEK_SET);

		mReferenceSequences.resize(mNumRefSeqs);
		mRefSeqLUT = new char*[mNumRefSeqs];

		unsigned int currentRefSeq = 0;
		vector<ReferenceSequence>::iterator rsIter;
		for(rsIter = mReferenceSequences.begin(); rsIter != mReferenceSequences.end(); ++rsIter, ++currentRefSeq) {

			// REFERENCE_SEQ_NAME_LEN[1]                0 -  0 
			// REFERENCE_SEQ_SPECIES_LEN[1]             1 -  1
			// REFERENCE_SEQ_GENOME_ASSEMBLY_ID_LEN[1]  2 -  2
			// REFERENCE_SEQ_URI_LEN[1]                 3 -  3
			// REFERENCE_SEQ_NUM_BASES[4]               4 -  7
			// REFERENCE_SEQ_SEQ_OFFSET[8]              8 - 15
			// REFERENCE_SEQ_MD5[16]                   16 - 31
			// REFERENCE_SEQ_NAME[X]                   32 - XX
			// REFERENCE_SEQ_SPECIES[X]
			// REFERENCE_SEQ_GENOME_ASSEMBLY_ID[X]
			// REFERENCE_SEQ_URI[X]

			// read the name length
			const unsigned char nameLen = fgetc(mInStream);

			// read the species length
			const unsigned char speciesLen = fgetc(mInStream);

			// read the genome assembly id length
			const unsigned char genomeAssemblyIDLen = fgetc(mInStream);

			// read the uri length
			const unsigned char uriLen = fgetc(mInStream);

			// read the number of bases
			fread((char*)&rsIter->NumBases, SIZEOF_INT, 1, mInStream);

			// write the number of aligned reads
			fread((char*)&rsIter->NumAligned, SIZEOF_UINT64, 1, mInStream);

			// read the MD5 checksum
			rsIter->MD5.resize(32);
			char* pBuffer = (char*)rsIter->MD5.data();
			fread(pBuffer, 32, 1, mInStream);

			// read the reference name
			rsIter->Name.resize(nameLen);
			pBuffer = (char*)rsIter->Name.data();
			fread(pBuffer, nameLen, 1, mInStream);

			mRefSeqLUT[currentRefSeq] = new char[nameLen + 1];
			memcpy(mRefSeqLUT[currentRefSeq], pBuffer, nameLen);
			mRefSeqLUT[currentRefSeq][nameLen] = 0;

			// read the species name
			if(speciesLen > 0) {
				rsIter->Species.resize(speciesLen);
				pBuffer = (char*)rsIter->Species.data();
				fread(pBuffer, speciesLen, 1, mInStream);
			}

			// read the genome assembly ID
			if(genomeAssemblyIDLen > 0) {
				rsIter->GenomeAssemblyID.resize(genomeAssemblyIDLen);
				pBuffer = (char*)rsIter->GenomeAssemblyID.data();
				fread(pBuffer, genomeAssemblyIDLen, 1, mInStream);
			}

			// read the URI
			if(uriLen > 0) {
				rsIter->URI.resize(uriLen);
				pBuffer = (char*)rsIter->URI.data();
				fread(pBuffer, uriLen, 1, mInStream);
			}

			// retrieve the number of reference sequence tags
			const unsigned char numReferenceSequenceTags = (unsigned char)fgetc(mInStream);

			if(numReferenceSequenceTags != 0) {
				printf("ERROR: Found reference sequence tags, but support for reference sequence tags has not been implemented yet.\n");
				exit(1);
			}

			//// DEBUG
			//cout << "# bases:                " << rsIter->NumBases << endl;
			//cout << "md5:                    " << rsIter->MD5 << endl;
			//cout << "name:                   " << rsIter->Name << endl;
			//cout << "species:                " << rsIter->Species << endl;
			//cout << "genome assembly ID:     " << rsIter->GenomeAssemblyID << endl;
			//cout << "URI:                    " << rsIter->URI << endl;
		}

		// ================================
		// read the reference sequence gaps
		// ================================

		//CFastLZIO fio;
		//if(mReferenceGapOffset != 0) {

		//	// jump to the reference gap location
		//	fseek64(mInStream, mReferenceGapOffset, SEEK_SET);

		//	// read the reference gaps vector
		//	fio.Read(mBuffer, mBufferLen, mInStream);

		//	unsigned int bufferOffset = 0;
		//	vector<GapInfo>::iterator gvIter;
		//	vector<vector<GapInfo> >::iterator rsgIter;
		//	for(rsgIter = mRefSeqGaps.begin(); rsgIter != mRefSeqGaps.end(); ++rsgIter) {

		//		// retrieve the number of gaps for this reference sequence
		//		unsigned int numGaps = 0;
		//		memcpy((char*)&numGaps, mBuffer + bufferOffset, SIZEOF_INT);
		//		bufferOffset += SIZEOF_INT;

		//		// pre-allocate the reference gap vector
		//		rsgIter->resize(numGaps);

		//		for(gvIter = rsgIter->begin(); gvIter != rsgIter->end(); ++gvIter) {

		//			// retrieve the reference gap position
		//			memcpy((char*)&gvIter->Position, mBuffer + bufferOffset, SIZEOF_INT);
		//			bufferOffset += SIZEOF_INT;

		//			// retrieve the reference gap length
		//			memcpy((char*)&gvIter->Length, mBuffer + bufferOffset, SIZEOF_SHORT);
		//			bufferOffset += SIZEOF_SHORT;
		//		}
		//	}
		//}

		// restore our file position
		Rewind();
	}

	// deserializes each alignment and stores them in the supplied vector
	void CAlignmentReader::ReadAlignments(vector<Alignment>& alignments, const bool isLongRead, const bool isPairedInSequencing, const bool isResolvedAsPair, const unsigned int readGroupCode) {
		vector<Alignment>::iterator alIter;
		for(alIter = alignments.begin(); alIter != alignments.end(); ++alIter) {
			ReadAlignment(*alIter, isLongRead, isPairedInSequencing, isResolvedAsPair);
			alIter->ReadGroupCode = readGroupCode;
		}
	}

	// deserialize the alignment
	void CAlignmentReader::ReadAlignment(Alignment& al, const bool isLongRead, const bool isPairedInSequencing, const bool isResolvedAsPair) {

		// get the reference sequence start position
		memcpy((char*)&al.ReferenceBegin, mBufferPtr, SIZEOF_INT);
		mBufferPtr += SIZEOF_INT;

		// get the reference sequence end position
		memcpy((char*)&al.ReferenceEnd, mBufferPtr, SIZEOF_INT);
		mBufferPtr += SIZEOF_INT;

		// get the reference sequence index
		memcpy((char*)&al.ReferenceIndex, mBufferPtr, SIZEOF_INT);
		mBufferPtr += SIZEOF_INT;
		al.ReferenceName = mRefSeqLUT[al.ReferenceIndex];

		// get the alignment quality
		al.Quality = (unsigned char)*mBufferPtr;
		++mBufferPtr;

		// get the alignment status flag
		const unsigned char status = (unsigned char)*mBufferPtr;
		++mBufferPtr;

		al.IsFirstMate             = false;
		al.IsReverseStrand     = false;
		al.IsMateReverseStrand = false;

		al.IsPairedEnd      = isPairedInSequencing;
		al.IsResolvedAsPair = isResolvedAsPair;

		if(isPairedInSequencing) {
			if((status & AF_IS_FIRST_MATE) != 0) al.IsFirstMate = true;
			if(isResolvedAsPair && ((status & AF_IS_MATE_REVERSE_STRAND) != 0)) al.IsMateReverseStrand = true;
		}

		if((status & AF_IS_REVERSE_STRAND) != 0) al.IsReverseStrand = true;
		if((status & AF_WAS_RESCUED)       != 0) al.WasRescued          = true;

		// get the number of mismatches
		memcpy((char*)&al.NumMismatches, mBufferPtr, SIZEOF_SHORT);
		mBufferPtr += SIZEOF_SHORT;

		// get mate pair information
		if(isResolvedAsPair) {

			// get the mate reference sequence start position
			memcpy((char*)&al.MateReferenceBegin, mBufferPtr, SIZEOF_INT);
			mBufferPtr += SIZEOF_INT;

			// get the mate reference sequence end position
			memcpy((char*)&al.MateReferenceEnd, mBufferPtr, SIZEOF_INT);
			mBufferPtr += SIZEOF_INT;

			// get the mate reference sequence index
			memcpy((char*)&al.MateReferenceIndex, mBufferPtr, SIZEOF_INT);
			mBufferPtr += SIZEOF_INT;

		} else {

			al.MateReferenceBegin = 0;
			al.MateReferenceEnd   = 0;
			al.MateReferenceIndex = 0;
		}

		unsigned short pairwiseLength = 0;

		if(isLongRead) {

			// get the pairwise length
			memcpy((char*)&pairwiseLength, mBufferPtr, SIZEOF_SHORT);
			mBufferPtr += SIZEOF_SHORT;

			// get the query begin
			memcpy((char*)&al.QueryBegin, mBufferPtr, SIZEOF_SHORT);
			mBufferPtr += SIZEOF_SHORT;

			// get the query end
			memcpy((char*)&al.QueryEnd, mBufferPtr, SIZEOF_SHORT);
			mBufferPtr += SIZEOF_SHORT;

		} else {

			// get the pairwise length
			pairwiseLength = (unsigned char)*mBufferPtr;
			++mBufferPtr;

			// get the query begin
			al.QueryBegin = (unsigned char)*mBufferPtr;
			++mBufferPtr;

			// get the query end
			al.QueryEnd = (unsigned char)*mBufferPtr;
			++mBufferPtr;
		}

		// retrieve the packed pairwise alignment
		string packed;
		packed.resize(pairwiseLength);
		memcpy((void*)packed.data(), mBufferPtr, pairwiseLength);
		mBufferPtr += pairwiseLength;

		// unpack the pairwise query bases
		CSequenceUtilities::Unpack(packed, al.Reference, al.Query);

		// get the pairwise query base qualities
		const unsigned short bqLength = al.QueryEnd - al.QueryBegin + 1;
		al.BaseQualities.resize(bqLength);
		memcpy((void*)al.BaseQualities.data(), mBufferPtr, bqLength);
		mBufferPtr += bqLength;

		// read the number of tags present in this alignment
		const unsigned char numTags = (unsigned char)*mBufferPtr;
		++mBufferPtr;

		if(numTags != 0) {
			cout << "ERROR: Tags have not been implemented yet." << endl;
			exit(1);
		}

		// DEBUG
		//cout << "reference index: " << al.ReferenceIndex << ", begin: " << al.ReferenceBegin << ", end: " << al.ReferenceEnd << endl;
		//cout << "query begin: " << al.QueryBegin << ", end: " << al.QueryEnd << ", pairwise length: " << pairwiseLength << endl << endl;
	}

	// reads a new compressed partition (returns false if EOF occurs)
	bool CAlignmentReader::ReadPartition(void) {

		// read the uncompressed partition entry size
		unsigned int uncompressedSize = 0;
		fread((char*)&uncompressedSize, SIZEOF_INT, 1, mInStream);

		if(feof(mInStream)) return false;

		// read the compressed partition entry size
		int compressedSize = 0;
		fread((char*)&compressedSize, SIZEOF_INT, 1, mInStream);

		// read the partition member size
		mPartitionMembers = 0;
		fread((char*)&mPartitionSize, SIZEOF_SHORT, 1, mInStream);

		// check the compression buffer size
		CMemoryUtilities::CheckBufferSize(mCompressionBuffer, mCompressionBufferLen, compressedSize);
		CMemoryUtilities::CheckBufferSize(mBuffer, mBufferLen, uncompressedSize);

		// read and uncompress the partition
		int numBytesRead = fread(mCompressionBuffer, 1, compressedSize, mInStream);

		if(numBytesRead != compressedSize) {
			cout << "ERROR: Tried to read " << compressedSize << " bytes, but received only " << numBytesRead << " bytes (" << mInputFilename << ") [read the partition: LoadNextRead]" << endl;
			exit(1);
		}

		int result = fastlz_decompress(mCompressionBuffer, compressedSize, mBuffer, mBufferLen);

		if(result == 0) {
			cout << "ERROR: Unable to properly uncompress the current data partition." << endl;
			exit(1);
		}

		// set the buffer pointer
		mBufferPtr = mBuffer;

		return true;
	}

	// reads the tag from disk
	void CAlignmentReader::ReadTag(Tag& tag) {

		// read the tag ID
		tag.ID = fgetc(mInStream);

		// read the tag type
		TagType tagType = fgetc(mInStream);
		tag.Type = tagType;

		// read the data
		unsigned short stringLength;

		switch(tagType) {
			case TT_CHAR:
				tag.Char = fgetc(mInStream);
				break;
			case TT_DOUBLE:
				fread((char*)&tag.Double, SIZEOF_DOUBLE, 1, mInStream);
				break;
			case TT_FLOAT:
				fread((char*)&tag.Float, SIZEOF_FLOAT, 1, mInStream);
				break;
			case TT_INT16:
				fread((char*)&tag.Int16, SIZEOF_SHORT, 1, mInStream);
				break;
			case TT_INT32:
				fread((char*)&tag.Int32, SIZEOF_INT, 1, mInStream);
				break;
			case TT_INT64:
				fread((char*)&tag.Int64, SIZEOF_UINT64, 1, mInStream);
				break;
			case TT_STRING:
				fread((char*)&stringLength, SIZEOF_SHORT, 1, mInStream);
				fread(tag.String, stringLength, 1, mInStream);
				break;
			case TT_UCHAR:
				tag.UChar = fgetc(mInStream);
				break;
			case TT_UINT16:
				fread((char*)&tag.UInt16, SIZEOF_SHORT, 1, mInStream);
				break;
			case TT_UINT32:
				fread((char*)&tag.UInt32, SIZEOF_INT, 1, mInStream);
				break;
			case TT_UINT64:
				fread((char*)&tag.UInt64, SIZEOF_UINT64, 1, mInStream);
				break;
			default:
				printf("ERROR: Unknown tag storage type found: %c [%u]\n", tagType, tagType);
				exit(1);
		}
	}

	// sets the file pointer to the beginning of the read data
	void CAlignmentReader::Rewind(void) {
		fseek64(mInStream, mReadsOffset, SEEK_SET);
		mCurrentRead      = 0;
		mPartitionMembers = 0;
		mPartitionSize    = 0;
	}

	// ================
	// CAlignmentWriter
	// ================

	// constructor
	CAlignmentWriter::CAlignmentWriter(void)
		: mIsOpen(false)
		, mOutStream(NULL)
		, mNumReads(0)
		, mNumBases(0)
		, mBuffer(NULL)
		, mBufferLen(10485760)
		, mBufferPosition(0)
		, mCompressionBuffer(NULL)
		, mCompressionBufferLen(0)
		, mPartitionSize(20000)
		, mPartitionMembers(0)
		, mpRefGapVector(NULL)
		, mStatus(AS_UNKNOWN)
		, mIsPairedEndArchive(false)
		, mLastReferenceIndex(0)
		, mLastReferencePosition(0)
		, mStoreIndex(false)
	{
		// set the buffer threshold
		mBufferThreshold = mBufferLen - MEMORY_BUFFER_SIZE;

		// initialize the read and index buffer
		try {
			mBuffer = new unsigned char[mBufferLen];
		} catch(bad_alloc) {
			cout << "ERROR: Unable to allocate enough memory for the alignment output buffer." << endl;
			exit(1);
		}
	}

	// destructor
	CAlignmentWriter::~CAlignmentWriter(void) {
		if(mIsOpen)            Close();
		if(mBuffer)            delete [] mBuffer;
		if(mCompressionBuffer) delete [] mCompressionBuffer;
	}

	// adds a header tag
	void CAlignmentWriter::AddHeaderTag(const Tag& tag) {

		// we can only add header tags before opening the file
		if(mIsOpen) {
			printf("ERROR: Header tags can only be added before opening the alignment archive.\n");
			exit(1);
		}

		map<unsigned char, Tag>::iterator htIter = mHeaderTags.find(tag.ID);
		if(htIter == mHeaderTags.end()) {
			mHeaderTags[tag.ID] = tag;
		} else htIter->second = tag;
	}

	// checks the buffer
	void CAlignmentWriter::AdjustBuffer(void) {

		// allocate a new buffer
		unsigned int newBufferLen = mBufferLen << 1;
		unsigned char* newBuffer = NULL;

		try {
			newBuffer = new unsigned char[newBufferLen];
		} catch(bad_alloc) {
			cout << "ERROR: Unable to reallocate enough memory for the alignment output buffer." << endl;
			exit(1);
		}

		// copy the old data and destroy the old buffer
		memcpy(newBuffer, mBuffer, mBufferLen);
		delete [] mBuffer;

		// repoint the new buffer
		mBuffer          = newBuffer;
		mBufferLen       = newBufferLen;
		mBufferThreshold = newBufferLen - MEMORY_BUFFER_SIZE;
	}

	// closes the alignment archive
	void CAlignmentWriter::Close(void) {

		// prevent the archive from being updated elsewhere
		mIsOpen = false;

		// flush the buffer
		if(mPartitionMembers > 0) WritePartition();

		// =======================================
		// save the reference sequence information
		// =======================================

		const off_type referenceOffset = ftell64(mOutStream);

		// write the reference sequence dictionary
		vector<ReferenceSequence>::const_iterator rsIter;
		for(rsIter = mReferenceSequences.begin(); rsIter != mReferenceSequences.end(); rsIter++) {

			// REFERENCE_SEQ_NAME_LEN[1]                0 -  0 
			// REFERENCE_SEQ_SPECIES_LEN[1]             1 -  1
			// REFERENCE_SEQ_GENOME_ASSEMBLY_ID_LEN[1]  2 -  2
			// REFERENCE_SEQ_URI_LEN[1]                 3 -  3
			// REFERENCE_SEQ_NUM_BASES[4]               4 -  7
			// REFERENCE_SEQ_SEQ_OFFSET[8]              8 - 15
			// REFERENCE_SEQ_MD5[16]                   16 - 31
			// REFERENCE_SEQ_NAME[X]                   32 - XX
			// REFERENCE_SEQ_SPECIES[X]
			// REFERENCE_SEQ_GENOME_ASSEMBLY_ID[X]
			// REFERENCE_SEQ_URI[X]

			// get the string lengths
			unsigned int nameLen             = (unsigned int)rsIter->Name.size();
			unsigned int speciesLen          = (unsigned int)rsIter->Species.size();
			unsigned int genomeAssemblyIDLen = (unsigned int)rsIter->GenomeAssemblyID.size();
			unsigned int uriLen              = (unsigned int)rsIter->URI.size();

			// write the name length
			fputc((unsigned char)nameLen, mOutStream);

			// write the species length
			fputc((unsigned char)speciesLen, mOutStream);

			// write the genome assembly id length
			fputc((unsigned char)genomeAssemblyIDLen, mOutStream);

			// write the URI length
			fputc((unsigned char)uriLen, mOutStream);

			// write the number of bases
			fwrite((char*)&rsIter->NumBases, SIZEOF_INT, 1, mOutStream);

			// write the number of aligned reads
			fwrite((char*)&rsIter->NumAligned, SIZEOF_UINT64, 1, mOutStream);

			// write the MD5 checksum
			fwrite(rsIter->MD5.data(), 32, 1, mOutStream);

			// write the reference name
			fwrite(rsIter->Name.data(), nameLen, 1, mOutStream);

			// write the species name
			if(speciesLen > 0) fwrite(rsIter->Species.data(), speciesLen, 1, mOutStream);

			// write the genome assembly ID
			if(genomeAssemblyIDLen > 0) fwrite(rsIter->GenomeAssemblyID.data(), genomeAssemblyIDLen, 1, mOutStream);

			// write the URI
			if(uriLen > 0) fwrite(rsIter->URI.data(), uriLen, 1, mOutStream);

			// write the number of reference sequence tags (hard coded as 0 for now)
			fputc(0, mOutStream);
		}

		// ==============
		// save the index
		// ==============

		CFastLZIO fio;
		off_type indexFileOffset = 0;
		if(mStoreIndex) {

			// get the current file offset
			indexFileOffset = ftell64(mOutStream);

			// calculate how much buffer space we'll need
			const unsigned int numIndexEntries = mIndex.size();
			const unsigned int requestedBytes = numIndexEntries * (SIZEOF_UINT64 + SIZEOF_INT + SIZEOF_SHORT) + SIZEOF_INT;
			CMemoryUtilities::CheckBufferSize(mBuffer, mBufferLen, requestedBytes);

			// store the number of index entries
			fwrite((char*)&numIndexEntries, SIZEOF_INT, 1, mOutStream);

			unsigned int bufferOffset = 0;
			vector<IndexEntry>::const_iterator indexIter;
			for(indexIter = mIndex.begin(); indexIter != mIndex.end(); indexIter++) {

				// store the reference index
				memcpy(mBuffer + bufferOffset, (char*)&indexIter->ReferenceIndex, SIZEOF_SHORT);
				bufferOffset += SIZEOF_SHORT;

				// store the reference position
				memcpy(mBuffer + bufferOffset, (char*)&indexIter->Position, SIZEOF_INT);
				bufferOffset += SIZEOF_INT;

				// store the file offset
				memcpy(mBuffer + bufferOffset, (char*)&indexIter->Offset, SIZEOF_UINT64);
				bufferOffset += SIZEOF_UINT64;
			}

			// write the index
			fio.Write((char*)mBuffer, bufferOffset, mOutStream);
		}

		// =================
		// update the header
		// =================

		// update the number of reads in the archive
		fseek64(mOutStream, NUM_READS_OFFSET, SEEK_SET);
		fwrite((char*)&mNumReads, SIZEOF_UINT64, 1, mOutStream);

		// update the number of bases in the archive
		fwrite((char*)&mNumBases, SIZEOF_UINT64, 1, mOutStream);

		// update the references offset in the archive
		fwrite((char*)&referenceOffset, SIZEOF_UINT64, 1, mOutStream);

		// update the reference gap offset in the archive
		off_type referenceGapFileOffset = 0;
		fwrite((char*)&referenceGapFileOffset, SIZEOF_UINT64, 1, mOutStream);

		// update the index offset in the archive		
		fwrite((char*)&indexFileOffset, SIZEOF_UINT64, 1, mOutStream);

		// DEBUG
		//cout << endl;
		//cout << "mNumReads:              " << mNumReads << endl;
		//cout << "mNumBases:              " << mNumBases << endl;
		//cout << "referenceOffset:        " << referenceOffset << endl;
		//cout << "referenceGapFileOffset: " << referenceGapFileOffset << endl;
		//cout << "indexFileOffset:        " << indexFileOffset << endl << endl;

		// close the file stream
		fclose(mOutStream);
	}

	// retrieves the number of bases written
	uint64_t CAlignmentWriter::GetNumBases(void) const {
		return mNumBases;
	}

	// retrieves the number of reads written
	uint64_t CAlignmentWriter::GetNumReads(void) const {
		return mNumReads;
	}

	// opens the alignment archive
	void CAlignmentWriter::Open(const string& filename, const vector<ReferenceSequence>& referenceSequences, const vector<ReadGroup>& readGroups, const AlignmentStatus as) {

		if(mIsOpen) {
			cout << "ERROR: An attempt was made to open an already open alignment archive." << endl;
			exit(1);
		}

		mOutputFilename = filename;

		mOutStream = fopen(filename.c_str(), "wb");
		if(!mOutStream) {
			cout << "ERROR: Could not open the compressed alignment archive (" << mOutputFilename << ") for writing." << endl;
			exit(1);
		}

		mIsOpen = true;

		// initialization
		mBufferPosition   = 0;
		mPartitionMembers = 0;
		mStatus           = as;

		// copy the reference sequence statistics
		// N.B. we copy these because g++ was making shallow copies before, this is only a temporary fix
		mNumRefSeqs = (unsigned int)referenceSequences.size();
		mReferenceSequences.resize(mNumRefSeqs);

		vector<ReferenceSequence>::const_iterator rsIter = referenceSequences.begin();
		vector<ReferenceSequence>::iterator crsIter;
		for(crsIter = mReferenceSequences.begin(); crsIter != mReferenceSequences.end(); crsIter++, rsIter++) {
			crsIter->GenomeAssemblyID = rsIter->GenomeAssemblyID;
			crsIter->MD5              = rsIter->MD5;
			crsIter->Name             = rsIter->Name;
			crsIter->NumAligned       = 0;
			crsIter->NumBases         = rsIter->NumBases;
			crsIter->Species          = rsIter->Species;
			crsIter->URI              = rsIter->URI;
		}

		// create our composite sequencing technology
		unsigned int numReadGroups = (unsigned int)readGroups.size();
		SequencingTechnologies st = ST_UNKNOWN;

		vector<ReadGroup>::const_iterator rgIter;
		for(rgIter = readGroups.begin(); rgIter != readGroups.end(); rgIter++) st |= rgIter->SequencingTechnology;

		// ================
		// write the header
		// ================

		// MOSAIK_SIGNATURE[6]	   0  -  5
		// STATUS[1]               6  -  6
		// SEQUENCE_TECHNOLOGY[2]  7  -  8
		// ARCHIVE_DATE[8]		   9  - 16
		// NUM_REFERENCE_SEQS[4]   17 - 20
		// NUM_READ_GROUPS[4]      21 - 24
		// NUM_READS[8]            25 - 32
		// NUM_BASES[8]            33 - 40
		// REFERENCES_OFFSET[8]    41 - 48
		// REFERENCE_GAP_OFFSET[8] 49 - 57
		// INDEX_OFFSET[8]         58 - 63
		// NUM_READ_GROUP_TAGS[1]  64 - 64
		// READ_GROUPS[*]

		// NB: the following blocks occur at the end of the file
		// REFERENCE_SEQS[*]
		// REFERENCE_GAPS[*]
		// INDEX[*]

		// write the MOSAIK signature
		const unsigned char SIGNATURE_LENGTH = 6;
		const char* MOSAIK_SIGNATURE = "MSKAA\4";
		fwrite(MOSAIK_SIGNATURE, SIGNATURE_LENGTH, 1, mOutStream);

		// write the alignment status
		fputc((unsigned char)as, mOutStream);
		if((mStatus & AS_SORTED_ALIGNMENT) != 0) mStoreIndex         = true;
		if((mStatus & AS_PAIRED_END_READ)  != 0) mIsPairedEndArchive = true;

		// write the sequencing technology
		fwrite((char*)&st, SIZEOF_SHORT, 1, mOutStream);

		// write the archive date
		uint64_t currentTime = CTimeSupport::GetSystemTime();
		fwrite((char*)&currentTime, SIZEOF_UINT64, 1, mOutStream);

		// write the number of reference sequences
		fwrite((char*)&mNumRefSeqs, SIZEOF_INT, 1, mOutStream);

		// write the number of read groups
		fwrite((char*)&numReadGroups, SIZEOF_INT, 1, mOutStream);

		// skip the # of reads, # of bases, references offset, reference gap offset, index offset
		fseek64(mOutStream, 2 * SIZEOF_UINT64 + 3 * SIZEOF_OFF_TYPE, SEEK_CUR);

		// write the number of header tags (hard coded as 0 for now)
		const unsigned char numHeaderTags = mHeaderTags.size();
		fputc(numHeaderTags, mOutStream);

		// write the header tags
		if(numHeaderTags != 0) {
			map<unsigned char, Tag>::const_iterator htIter;
			for(htIter = mHeaderTags.begin(); htIter != mHeaderTags.end(); htIter++) {
				WriteTag(htIter);
			}
		}

		// write the read groups
		for(rgIter = readGroups.begin(); rgIter != readGroups.end(); rgIter++) {

			// write the metadata string lengths
			const unsigned char centerNameLen   = (unsigned char)rgIter->CenterName.size();
			const unsigned char libraryNameLen  = (unsigned char)rgIter->LibraryName.size();
			const unsigned char platformUnitLen = (unsigned char)rgIter->PlatformUnit.size();
			const unsigned char readGroupIDLen  = (unsigned char)rgIter->ReadGroupID.size();
			const unsigned char sampleNameLen   = (unsigned char)rgIter->SampleName.size();
			const unsigned short descriptionLen = (unsigned short)rgIter->Description.size();

			fputc(centerNameLen,   mOutStream);
			fputc(libraryNameLen,  mOutStream);
			fputc(platformUnitLen, mOutStream);
			fputc(readGroupIDLen,  mOutStream);
			fputc(sampleNameLen,   mOutStream);
			fwrite((char*)&descriptionLen, SIZEOF_SHORT, 1, mOutStream);
			fwrite((char*)&rgIter->SequencingTechnology, SIZEOF_SHORT, 1, mOutStream);
			fwrite((char*)&rgIter->MedianFragmentLength, SIZEOF_INT, 1, mOutStream);

			// write the metadata strings
			fwrite(rgIter->CenterName.c_str(),   centerNameLen,   1, mOutStream);
			fwrite(rgIter->Description.c_str(),  descriptionLen,  1, mOutStream);
			fwrite(rgIter->LibraryName.c_str(),  libraryNameLen,  1, mOutStream);
			fwrite(rgIter->PlatformUnit.c_str(), platformUnitLen, 1, mOutStream);
			fwrite(rgIter->ReadGroupID.c_str(),  readGroupIDLen,  1, mOutStream);
			fwrite(rgIter->SampleName.c_str(),   sampleNameLen,   1, mOutStream);

			// write the number of read group tags (hard coded as 0 for now)
			fputc(0, mOutStream);
		}
	}

	// saves the read to the alignment archive
	void CAlignmentWriter::SaveAlignedRead(const Mosaik::AlignedRead& ar) {

		// check the memory buffer
		if(mBufferPosition > mBufferThreshold) AdjustBuffer();

		// initialize
		const unsigned int numMate1Alignments = ar.Mate1Alignments.size();
		const unsigned int numMate2Alignments = ar.Mate2Alignments.size();

		const bool haveMate1 = (numMate1Alignments != 0 ? true : false);
		const bool haveMate2 = (numMate2Alignments != 0 ? true : false);

		// check if this is a long read
		const bool isLongRead = ar.IsLongRead;

		// derive our read status
		unsigned char readStatus = RF_UNKNOWN;

		if(haveMate1)      readStatus |= RF_HAVE_MATE1;
		if(haveMate2)      readStatus |= RF_HAVE_MATE2;
		if(isLongRead)     readStatus |= RF_IS_LONG_READ;
		if(ar.IsPairedEnd) readStatus |= RF_IS_PAIRED_IN_SEQUENCING;

		// write the read header
		WriteReadHeader(ar.Name, ar.ReadGroupCode, readStatus, numMate1Alignments, numMate2Alignments);

		// ===============================
		// serialize each mate 1 alignment
		// ===============================

		vector<Alignment>::const_iterator alIter;
		const Alignment *pAl = NULL, *pAlBegin = NULL;

		if(haveMate1) {
			pAlBegin = &ar.Mate1Alignments[0];
			for(alIter = ar.Mate1Alignments.begin(); alIter != ar.Mate1Alignments.end(); ++alIter) {
				pAl = pAlBegin + (alIter - ar.Mate1Alignments.begin());
				WriteAlignment(pAl, isLongRead, ar.IsPairedEnd, true, false);
			}
		}

		// ===============================
		// serialize each mate 2 alignment
		// ===============================

		if(haveMate2) {
			pAlBegin = &ar.Mate2Alignments[0];
			for(alIter = ar.Mate2Alignments.begin(); alIter != ar.Mate2Alignments.end(); ++alIter) {
				pAl = pAlBegin + (alIter - ar.Mate2Alignments.begin());
				WriteAlignment(pAl, isLongRead, ar.IsPairedEnd, false, false);
			}
		}

		// flush the buffer
		mPartitionMembers++;
		if(mPartitionMembers >= mPartitionSize) WritePartition();

		// increment the read counter
		mNumReads++;
	}

	// serializes the specified alignment
	void CAlignmentWriter::WriteAlignment(const Alignment* pAl, const bool isLongRead, const bool isPairedEnd, const bool isFirstMate, const bool isResolvedAsPair) {

		// check the memory buffer
		if(mBufferPosition > mBufferThreshold) AdjustBuffer();

		// store the reference sequence start position
		memcpy(mBuffer + mBufferPosition, (char*)&pAl->ReferenceBegin, SIZEOF_INT);
		mBufferPosition += SIZEOF_INT;

		// store the reference sequence end position
		memcpy(mBuffer + mBufferPosition, (char*)&pAl->ReferenceEnd, SIZEOF_INT);
		mBufferPosition += SIZEOF_INT;

		// store the reference sequence index
		mReferenceSequences[pAl->ReferenceIndex].NumAligned++;
		memcpy(mBuffer + mBufferPosition, (char*)&pAl->ReferenceIndex, SIZEOF_INT);
		mBufferPosition += SIZEOF_INT;

		// store the alignment quality
		mBuffer[mBufferPosition++] = pAl->Quality;

		// store the alignment status flag
		unsigned char status = AF_UNKNOWN;

		if(isPairedEnd && isFirstMate)                       status |= AF_IS_FIRST_MATE;
		if(isPairedEnd && !isFirstMate)                      status |= AF_IS_SECOND_MATE;
		if(pAl->IsReverseStrand)                         status |= AF_IS_REVERSE_STRAND;
		if(isResolvedAsPair && pAl->IsMateReverseStrand) status |= AF_IS_MATE_REVERSE_STRAND;
		if(pAl->WasRescued)                                  status |= AF_WAS_RESCUED;

		// not really sure how this applies to single-end and paired-end reads. Disabling the flag for now.
		//if(!isPrimaryAlignment)                              status |= AF_IS_NOT_PRIMARY;

		mBuffer[mBufferPosition++] = status;

		// store the number of mismatches
		memcpy(mBuffer + mBufferPosition, (char*)&pAl->NumMismatches, SIZEOF_SHORT);
		mBufferPosition += SIZEOF_SHORT;

		// get the pairwise length
		const unsigned short pairwiseLength = (unsigned short)pAl->Reference.size();

		// add mate pair information
		if(isResolvedAsPair) {

			// store the mate reference sequence start position
			memcpy(mBuffer + mBufferPosition, (char*)&pAl->MateReferenceBegin, SIZEOF_INT);
			mBufferPosition += SIZEOF_INT;

			// store the mate reference sequence end position
			memcpy(mBuffer + mBufferPosition, (char*)&pAl->MateReferenceEnd, SIZEOF_INT);
			mBufferPosition += SIZEOF_INT;

			// store the mate reference sequence index
			memcpy(mBuffer + mBufferPosition, (char*)&pAl->MateReferenceIndex, SIZEOF_INT);
			mBufferPosition += SIZEOF_INT;
		}

		if(isLongRead) {

			// store the pairwise length
			memcpy(mBuffer + mBufferPosition, (char*)&pairwiseLength, SIZEOF_SHORT);
			mBufferPosition += SIZEOF_SHORT;

			// store the query begin
			memcpy(mBuffer + mBufferPosition, (char*)&pAl->QueryBegin, SIZEOF_SHORT);
			mBufferPosition += SIZEOF_SHORT;

			// store the query end
			memcpy(mBuffer + mBufferPosition, (char*)&pAl->QueryEnd, SIZEOF_SHORT);
			mBufferPosition += SIZEOF_SHORT;

		} else {

			// store the pairwise length
			mBuffer[mBufferPosition++] = (unsigned char)pairwiseLength;

			// store the query begin
			mBuffer[mBufferPosition++] = (unsigned char)pAl->QueryBegin;

			// store the query end
			mBuffer[mBufferPosition++] = (unsigned char)pAl->QueryEnd;
		}

		// pack the pairwise-alignment string
		string packString;
		CSequenceUtilities::Pack(packString, pAl->Reference, pAl->Query);

		// store the packed pairwise alignment
		memcpy(mBuffer + mBufferPosition, packString.c_str(), pairwiseLength);
		mBufferPosition += pairwiseLength;

		// DEBUG
		//printf("\nReference: %s\n", pAl->Reference.c_str());
		//printf("Query:     %s\n", pAl->Query.c_str());
		////printf("Packed data: ");
		////for(unsigned short k = 0; k < pairwiseLength; ++k) {
		////	printf("%02x ", packString[k]);
		////}
		//const char* pRef = pAl->Reference.c_str();
		//const char* pQry = pAl->Query.c_str();
		//for(unsigned short k = 0; k < pairwiseLength; ++k) {
		//	if(pRef[k] != pQry[k]) {
		//		printf("ref: %c, query: %c, pack string: %02x\n", pRef[k], pQry[k], packString[k]);
		//	}
		//	//printf("%02x ", packString[k]);
		//}
		//printf("\n\n");
		////exit(1);

		// store the pairwise query base qualities
		const unsigned int bqLength = pAl->BaseQualities.size();
		memcpy(mBuffer + mBufferPosition, pAl->BaseQualities.c_str(), bqLength);
		mBufferPosition += bqLength;

		// write the number of tags present in this alignment (hard coded as 0 for now)
		mBuffer[mBufferPosition++] = 0;

		// update our statistics
		mNumBases += bqLength;

		// check the buffer
		if(mBufferPosition >= mBufferLen) {
			cout << endl << "ERROR: Buffer overrun detected when saving read. Used " << mBufferPosition << " bytes, but allocated " << mBufferLen << " bytes." << endl;
			exit(1);
		}
	}

	// write the read header to disk
	void CAlignmentWriter::WriteReadHeader(const string& readName, const unsigned int readGroupCode, const unsigned char readStatus, const unsigned int numMate1Alignments, const unsigned int numMate2Alignments) {

		// store the read name
		const unsigned char readNameLen = (unsigned char)readName.size();
		mBuffer[mBufferPosition++] = readNameLen;
		memcpy(mBuffer + mBufferPosition, readName.c_str(), readNameLen);
		mBufferPosition += readNameLen;

		// store the read group code
		memcpy(mBuffer + mBufferPosition, (char*)&readGroupCode, SIZEOF_INT);
		mBufferPosition += SIZEOF_INT;

		// store the read status flag
		mBuffer[mBufferPosition++] = readStatus;

		// store the number of mate 1 alignments
		if(numMate1Alignments != 0) {
			memcpy(mBuffer + mBufferPosition, (char*)&numMate1Alignments, SIZEOF_INT);
			mBufferPosition += SIZEOF_INT;
		}

		// store the number of mate 2 alignments
		if(numMate2Alignments != 0) {
			memcpy(mBuffer + mBufferPosition, (char*)&numMate2Alignments, SIZEOF_INT);
			mBufferPosition += SIZEOF_INT;
		}
	}

	// write partition to disk
	void CAlignmentWriter::WritePartition(void) {

		// store the partition index entry
		if(mStoreIndex) {
			IndexEntry ie;
			ie.Offset         = ftell64(mOutStream);
			ie.ReferenceIndex = mLastReferenceIndex;
			ie.Position       = mLastReferencePosition;
			mIndex.push_back(ie);
		}

		// check the compression buffer size
		unsigned int requestedSize = (unsigned int)(mBufferPosition * 1.05);
		CMemoryUtilities::CheckBufferSize(mCompressionBuffer, mCompressionBufferLen, requestedSize);

		// compress the partition
		int compressedSize = fastlz_compress_level(FASTLZ_BETTER_COMPRESSION, mBuffer, mBufferPosition, mCompressionBuffer);

		// write the uncompressed partition entry size
		fwrite((char*)&mBufferPosition, SIZEOF_INT, 1, mOutStream);

		// write the compressed partition entry size
		fwrite((char*)&compressedSize, SIZEOF_INT, 1, mOutStream);

		// write the partition member size
		fwrite((char*)&mPartitionMembers, SIZEOF_SHORT, 1, mOutStream);

		// write the partition
		fwrite(mCompressionBuffer, compressedSize, 1, mOutStream);

		mPartitionMembers = 0;
		mBufferPosition   = 0;
	}

	// writes the tag to disk
	void CAlignmentWriter::WriteTag(const map<unsigned char, Tag>::const_iterator& htIter) {

		// localize the tag type
		TagType tagType = htIter->second.Type;

		// write the tag ID
		fputc(htIter->second.ID, mOutStream);

		// write the tag type
		fputc(tagType, mOutStream);

		// write the data
		unsigned short stringLength;

		switch(tagType) {
			case TT_CHAR:
				fputc(htIter->second.Char, mOutStream);
				break;
			case TT_DOUBLE:
				fwrite((char*)&htIter->second.Double, SIZEOF_DOUBLE, 1, mOutStream);
				break;
			case TT_FLOAT:
				fwrite((char*)&htIter->second.Float, SIZEOF_FLOAT, 1, mOutStream);
				break;
			case TT_INT16:
				fwrite((char*)&htIter->second.Int16, SIZEOF_SHORT, 1, mOutStream);
				break;
			case TT_INT32:
				fwrite((char*)&htIter->second.Int32, SIZEOF_INT, 1, mOutStream);
				break;
			case TT_INT64:
				fwrite((char*)&htIter->second.Int64, SIZEOF_UINT64, 1, mOutStream);
				break;
			case TT_STRING:
				stringLength = (unsigned short)strlen(htIter->second.String);
				fwrite((char*)&stringLength, SIZEOF_SHORT, 1, mOutStream);
				fwrite(htIter->second.String, stringLength, 1, mOutStream);
				break;
			case TT_UCHAR:
				fputc(htIter->second.UChar, mOutStream);
				break;
			case TT_UINT16:
				fwrite((char*)&htIter->second.UInt16, SIZEOF_SHORT, 1, mOutStream);
				break;
			case TT_UINT32:
				fwrite((char*)&htIter->second.UInt32, SIZEOF_INT, 1, mOutStream);
				break;
			case TT_UINT64:
				fwrite((char*)&htIter->second.UInt64, SIZEOF_UINT64, 1, mOutStream);
				break;
			default:
				printf("ERROR: Unknown tag storage type found: %c\n", tagType);
				exit(1);
		}
	}

	// ================
	// CMemoryUtilities
	// ================

	// checks if the buffer is large enough to accommodate the requested size
	void CMemoryUtilities::CheckBufferSize(char* &pBuffer, unsigned int& bufferLen, const unsigned int requestedBytes) {
		try {
			if(requestedBytes > bufferLen) {
				bufferLen = requestedBytes + 10;
				delete [] pBuffer;
				pBuffer = new char[bufferLen];
			}
		} catch(bad_alloc) {
			cout << "ERROR: Out of memory when allocating " << requestedBytes << " bytes." << endl;
			exit(1);
		}
	}

	// checks if the buffer is large enough to accommodate the requested size
	void CMemoryUtilities::CheckBufferSize(unsigned char* &pBuffer, unsigned int& bufferLen, const unsigned int requestedBytes) {
		try {
			if(requestedBytes > bufferLen) {
				bufferLen = requestedBytes + 10;
				delete [] pBuffer;
				pBuffer = new unsigned char[bufferLen];
			}
		} catch(bad_alloc) {
			cout << "ERROR: Out of memory when allocating " << requestedBytes << " bytes." << endl;
			exit(1);
		}
	}

	// ============
	// CTimeSupport
	// ============

#ifdef WIN32
	const uint64_t CTimeSupport::SECS_BETWEEN_EPOCHS = 11644473600;
	const uint64_t CTimeSupport::SECS_TO_100NS       = 10000000;
#else
	const uint64_t CTimeSupport::SECS_BETWEEN_EPOCHS = 11644473600llu;
	const uint64_t CTimeSupport::SECS_TO_100NS       = 10000000llu;
#endif

	// returns the current time (UTC)
	uint64_t CTimeSupport::GetSystemTime(void) {

		uint64_t currentTime = 0;

#ifdef WIN32
		FILETIME ft;
		ULARGE_INTEGER ul;
		GetSystemTimeAsFileTime(&ft);
		ul.HighPart = ft.dwHighDateTime;
		ul.LowPart  = ft.dwLowDateTime;
		currentTime = ul.QuadPart;
#else
		currentTime = ConvertTimeT(time(NULL));
#endif

		return currentTime;
	}

	// converts a time_t variable to our 64-bit notation
	uint64_t CTimeSupport::ConvertTimeT(const time_t& timeT) {
		return (timeT + SECS_BETWEEN_EPOCHS) * SECS_TO_100NS;
	}

	// =========
	// CFastLZIO
	// =========

	// constructor
	CFastLZIO::CFastLZIO(void)
		: mBuffer(NULL)
		, mBufferLen(0)
	{
		try {
			mBufferLen = FASTLZ_IO_BUFFER_LEN;
			mBuffer = new char[mBufferLen];
		} catch(bad_alloc) {
			printf("ERROR: Unable to initialize the FastLZ buffer.\n");
			exit(1);
		}
	}

	// destructor
	CFastLZIO::~CFastLZIO(void) {
		Clear();
	}

	// clears the buffer
	void CFastLZIO::Clear(void) {
		mBufferLen = 0;
		if(mBuffer) {
			delete [] mBuffer;
			mBuffer = NULL;
		}
	}

	// our input method
	void CFastLZIO::Read(char* &buffer, unsigned int& bufferLen, FILE* stm) {

		// read the buffer length
		unsigned int newBufferLen;
		fread((char*)&newBufferLen, SIZEOF_INT, 1, stm);

		// allocate memory if necessary
		if(newBufferLen > bufferLen) {
			try {
				bufferLen = newBufferLen;
				if(buffer) delete [] buffer;
				buffer = new char[bufferLen + 1];
			} catch(bad_alloc) {
				printf("ERROR: Unable to initialize the uncompressed FastLZ buffer (Read).\n");
				exit(1);
			}
		}

		// calculate the number of blocks
		const unsigned short numBlocksRead = (unsigned short)ceil(bufferLen / (double)FASTLZ_IO_OUTPUT_BUFFER_LEN);

		// read each block
		char* pBuffer = buffer;
		int numCompressedBytes;

		for(unsigned int i = 0; i < numBlocksRead; ++i) {

			// read the block length
			fread((char*)&numCompressedBytes, SIZEOF_INT, 1, stm);

			// read the compressed block
			fread(mBuffer, numCompressedBytes, 1, stm);

			// uncompress the block
			int numUncompressedBytes = fastlz_decompress(mBuffer, numCompressedBytes, (void*)pBuffer, FASTLZ_IO_BUFFER_LEN);
			pBuffer += numUncompressedBytes;
		}

		// add the null termination
		*pBuffer = 0;
	}

	// our input method (STL string)
	void CFastLZIO::Read(string& s, FILE* stm) {

		// read the buffer length
		unsigned int bufferLen;
		fread((char*)&bufferLen, SIZEOF_INT, 1, stm);

		// handle an empty packet
		if(bufferLen == 0) {
			s.clear();
			return;
		}

		// calculate the number of blocks
		const unsigned short numBlocksRead = (unsigned short)ceil(bufferLen / (double)FASTLZ_IO_OUTPUT_BUFFER_LEN);

		// resize the string
		s.resize(bufferLen);

		// read each block
		char* pBuffer = (char*)s.data();
		int numCompressedBytes;

		for(unsigned int i = 0; i < numBlocksRead; ++i) {
			fread((char*)&numCompressedBytes, SIZEOF_INT, 1, stm);
			fread(mBuffer, numCompressedBytes, 1, stm);
			int numUncompressedBytes = fastlz_decompress(mBuffer, numCompressedBytes, (void*)pBuffer, FASTLZ_IO_BUFFER_LEN);
			pBuffer += numUncompressedBytes;
		}
	}

	// our output method
	void CFastLZIO::Write(const char* buffer, const unsigned int bufferLen, FILE* stm) {

		// write the buffer length
		fwrite((char*)&bufferLen, SIZEOF_INT, 1, stm);

		// calculate the number of blocks
		const unsigned short numBlocksWritten = (unsigned short)ceil(bufferLen / (double)FASTLZ_IO_OUTPUT_BUFFER_LEN);

		// write each block
		const char* pBuffer    = buffer;
		unsigned int bytesLeft = bufferLen;

		for(unsigned int i = 0; i < numBlocksWritten; ++i) {

			// compress the block
			unsigned int numUncompressedBytes = (bytesLeft > FASTLZ_IO_OUTPUT_BUFFER_LEN ? FASTLZ_IO_OUTPUT_BUFFER_LEN : bytesLeft);
			int numCompressedBytes = fastlz_compress_level(FASTLZ_BETTER_COMPRESSION, pBuffer, numUncompressedBytes, mBuffer);
			pBuffer   += numUncompressedBytes;
			bytesLeft -= numUncompressedBytes;

			// write the block length
			fwrite((char*)&numCompressedBytes, SIZEOF_INT, 1, stm);
			//printf("uncompressed bytes: %u, compressed bytes: %u\n", numUncompressedBytes, numCompressedBytes);

			// write the compressed block
			fwrite(mBuffer, numCompressedBytes, 1, stm);
		}
	}

	// =====
	// CSHA1
	// =====

	// generates the read group code based on the supplied text
	unsigned int CSHA1::GenerateReadGroupCode(const string& readGroupID, const string& sampleName) {

		// concatenate the two strings
		char concatenatedString[MAX_STRING_LENGTH];
		sprintf(concatenatedString, "%s%s", readGroupID.c_str(), sampleName.c_str());
		const unsigned int concatenatedStringLen = (unsigned int)strlen(concatenatedString);

		// initialize our SHA1 object
		CSHA1 sha1;

		// initialize our contexts
		SHA1_Data contexts[2];
		sha1.InitializeSHA1(contexts[0]);
		sha1.InitializeSHA1(contexts[1]);
		contexts[1].bytes = 1;

		// perform a SHA1 hash on the hash string
		sha1.UpdateSHA1(contexts[0], concatenatedString, concatenatedStringLen);
		sha1.UpdateSHA1(contexts[1], concatenatedString, concatenatedStringLen);

		// extract the final combined string from an array of hash private buffers
		unsigned int* pKey = (unsigned int*)sha1.mKey;
		sha1.FinalizeSHA1(contexts[0],     pKey, 5);
		sha1.FinalizeSHA1(contexts[1], pKey + 5, 3);

		//printf("SHA1: { ");
		//for(unsigned char m = 0; m < SHA_HASH_LENGTH; m++) printf("%02X ", sha1.mKey[m]);
		//printf("}\n");

		// copy the first 4 bytes from the hash
		unsigned int readGroupCode;
		memcpy((char*)&readGroupCode, sha1.mKey, 4);

		return readGroupCode;
	}

	// initialize the SHA values
	void CSHA1::InitializeSHA1(SHA1_Data& data) {

		// Set the h-vars to their initial values
		data.iv[0] = 0x67452301;
		data.iv[1] = 0xEFCDAB89;
		data.iv[2] = 0x98BADCFE;
		data.iv[3] = 0x10325476;
		data.iv[4] = 0xC3D2E1F0;

		// Initialise bit count
		data.bytes = 0;

		memset((char*)&data.key, 0, PGP_SHA_BLOCKBYTES);
	}

	// perform the SHA transformation
	void CSHA1::TransformSHA1(unsigned int* block, unsigned int* key) {

		// Set up first buffer
		register unsigned int A = block[0];
		register unsigned int B = block[1];
		register unsigned int C = block[2];
		register unsigned int D = block[3];
		register unsigned int E = block[4];
		register unsigned int t;

		// Heavy mangling, in 4 sub-rounds of 20 interations each.
		subRound( A, B, C, D, E, f1, K2, key[ 0] );
		subRound( E, A, B, C, D, f1, K2, key[ 1] );
		subRound( D, E, A, B, C, f1, K2, key[ 2] );
		subRound( C, D, E, A, B, f1, K2, key[ 3] );
		subRound( B, C, D, E, A, f1, K2, key[ 4] );
		subRound( A, B, C, D, E, f1, K2, key[ 5] );
		subRound( E, A, B, C, D, f1, K2, key[ 6] );
		subRound( D, E, A, B, C, f1, K2, key[ 7] );
		subRound( C, D, E, A, B, f1, K2, key[ 8] );
		subRound( B, C, D, E, A, f1, K2, key[ 9] );
		subRound( A, B, C, D, E, f1, K2, key[10] );
		subRound( E, A, B, C, D, f1, K2, key[11] );
		subRound( D, E, A, B, C, f1, K2, key[12] );
		subRound( C, D, E, A, B, f1, K2, key[13] );
		subRound( B, C, D, E, A, f1, K2, key[14] );
		subRound( A, B, C, D, E, f1, K2, key[15] );
		subRound( E, A, B, C, D, f1, K2, expand(key, 16) );
		subRound( D, E, A, B, C, f1, K2, expand(key, 17) );
		subRound( C, D, E, A, B, f1, K2, expand(key, 18) );
		subRound( B, C, D, E, A, f1, K2, expand(key, 19) );

		subRound( A, B, C, D, E, f2, K3, expand(key, 20) );
		subRound( E, A, B, C, D, f2, K3, expand(key, 21) );
		subRound( D, E, A, B, C, f2, K3, expand(key, 22) );
		subRound( C, D, E, A, B, f2, K3, expand(key, 23) );
		subRound( B, C, D, E, A, f2, K3, expand(key, 24) );
		subRound( A, B, C, D, E, f2, K3, expand(key, 25) );
		subRound( E, A, B, C, D, f2, K3, expand(key, 26) );
		subRound( D, E, A, B, C, f2, K3, expand(key, 27) );
		subRound( C, D, E, A, B, f2, K3, expand(key, 28) );
		subRound( B, C, D, E, A, f2, K3, expand(key, 29) );
		subRound( A, B, C, D, E, f2, K3, expand(key, 30) );
		subRound( E, A, B, C, D, f2, K3, expand(key, 31) );
		subRound( D, E, A, B, C, f2, K3, expand(key, 32) );
		subRound( C, D, E, A, B, f2, K3, expand(key, 33) );
		subRound( B, C, D, E, A, f2, K3, expand(key, 34) );
		subRound( A, B, C, D, E, f2, K3, expand(key, 35) );
		subRound( E, A, B, C, D, f2, K3, expand(key, 36) );
		subRound( D, E, A, B, C, f2, K3, expand(key, 37) );
		subRound( C, D, E, A, B, f2, K3, expand(key, 38) );
		subRound( B, C, D, E, A, f2, K3, expand(key, 39) );

		subRound( A, B, C, D, E, f3, K5, expand(key, 40) );
		subRound( E, A, B, C, D, f3, K5, expand(key, 41) );
		subRound( D, E, A, B, C, f3, K5, expand(key, 42) );
		subRound( C, D, E, A, B, f3, K5, expand(key, 43) );
		subRound( B, C, D, E, A, f3, K5, expand(key, 44) );
		subRound( A, B, C, D, E, f3, K5, expand(key, 45) );
		subRound( E, A, B, C, D, f3, K5, expand(key, 46) );
		subRound( D, E, A, B, C, f3, K5, expand(key, 47) );
		subRound( C, D, E, A, B, f3, K5, expand(key, 48) );
		subRound( B, C, D, E, A, f3, K5, expand(key, 49) );
		subRound( A, B, C, D, E, f3, K5, expand(key, 50) );
		subRound( E, A, B, C, D, f3, K5, expand(key, 51) );
		subRound( D, E, A, B, C, f3, K5, expand(key, 52) );
		subRound( C, D, E, A, B, f3, K5, expand(key, 53) );
		subRound( B, C, D, E, A, f3, K5, expand(key, 54) );
		subRound( A, B, C, D, E, f3, K5, expand(key, 55) );
		subRound( E, A, B, C, D, f3, K5, expand(key, 56) );
		subRound( D, E, A, B, C, f3, K5, expand(key, 57) );
		subRound( C, D, E, A, B, f3, K5, expand(key, 58) );
		subRound( B, C, D, E, A, f3, K5, expand(key, 59) );

		subRound( A, B, C, D, E, f4, K10, expand(key, 60) );
		subRound( E, A, B, C, D, f4, K10, expand(key, 61) );
		subRound( D, E, A, B, C, f4, K10, expand(key, 62) );
		subRound( C, D, E, A, B, f4, K10, expand(key, 63) );
		subRound( B, C, D, E, A, f4, K10, expand(key, 64) );
		subRound( A, B, C, D, E, f4, K10, expand(key, 65) );
		subRound( E, A, B, C, D, f4, K10, expand(key, 66) );
		subRound( D, E, A, B, C, f4, K10, expand(key, 67) );
		subRound( C, D, E, A, B, f4, K10, expand(key, 68) );
		subRound( B, C, D, E, A, f4, K10, expand(key, 69) );
		subRound( A, B, C, D, E, f4, K10, expand(key, 70) );
		subRound( E, A, B, C, D, f4, K10, expand(key, 71) );
		subRound( D, E, A, B, C, f4, K10, expand(key, 72) );
		subRound( C, D, E, A, B, f4, K10, expand(key, 73) );
		subRound( B, C, D, E, A, f4, K10, expand(key, 74) );
		subRound( A, B, C, D, E, f4, K10, expand(key, 75) );
		subRound( E, A, B, C, D, f4, K10, expand(key, 76) );
		subRound( D, E, A, B, C, f4, K10, expandx(key, 77) );
		subRound( C, D, E, A, B, f4, K10, expandx(key, 78) );
		subRound( B, C, D, E, A, f4, K10, expandx(key, 79) );

		// Build message digest
		block[0] += A;
		block[1] += B;
		block[2] += C;
		block[3] += D;
		block[4] += E;
	}

	// update SHA for a block of data
	void CSHA1::UpdateSHA1(SHA1_Data& data, void const *bufIn, size_t len) {

		unsigned char *buf = (unsigned char *)bufIn;

		// Update bitcount
		unsigned int i = (unsigned int)data.bytes % PGP_SHA_BLOCKBYTES;
		data.bytes += (unsigned int)len;

		// i is always less than PGP_SHA_BLOCKBYTES.
		if(PGP_SHA_BLOCKBYTES-i > len) {
			memmove((unsigned char *)data.key + i, buf, len);
			return;
		}

		if(i) {	// First chunk is an odd size
			memcpy((unsigned char *)data.key + i, buf, PGP_SHA_BLOCKBYTES - i);
			SwapBytes(data.key, (unsigned char *)data.key, NUM_SHA_KEY_WORDS);
			TransformSHA1(data.iv, data.key);
			buf += PGP_SHA_BLOCKBYTES-i;
			len -= PGP_SHA_BLOCKBYTES-i;
		}

		// Process data in 64-byte chunks
		while(len >= PGP_SHA_BLOCKBYTES) {
			SwapBytes(data.key, buf, NUM_SHA_KEY_WORDS);
			TransformSHA1(data.iv, data.key);
			buf += PGP_SHA_BLOCKBYTES;
			len -= PGP_SHA_BLOCKBYTES;
		}

		// Handle any remaining bytes of data.
		if(len) memmove(data.key, buf, len);
	}

	// final wrapup - pad to 64-byte boundary
	void CSHA1::FinalizeSHA1(SHA1_Data& data, unsigned int* key, const unsigned char numWords) {

		unsigned int i = (unsigned int)data.bytes % PGP_SHA_BLOCKBYTES;
		unsigned char *p = (unsigned char *)data.key + i; // First unused byte

		// Set the first char of padding to 0x80.  There is always room.
		*p++ = 0x80;

		// Bytes of padding needed to make 64 bytes (0..63)
		i = PGP_SHA_BLOCKBYTES - 1 - i;

		if(i < 8) {	// Padding forces an extra block
			memset(p, 0, i);
			SwapBytes(data.key, (unsigned char *)data.key, 16);
			TransformSHA1(data.iv, data.key);
			p = (unsigned char *)data.key;
			i = 64;
		}

		memset(p, 0, i - 8);
		SwapBytes(data.key, (unsigned char *)data.key, 14);

		// Append length in bits and transform
		data.key[14] = (unsigned int)(data.bytes >> 29);
		data.key[15] = (unsigned int)data.bytes << 3;
		TransformSHA1(data.iv, data.key);

		for(i = 0; i < NUM_SHA_IV_WORDS; i++) data.iv[i] = swap_int32(data.iv[i]);
		memcpy(key, data.iv, numWords * 4);
	}

	// Shuffle the bytes into big-endian order within words, as per the SHA spec.
	void CSHA1::SwapBytes(unsigned int* dest, unsigned char const *src, unsigned char words) {
		do {
			*dest++ = (unsigned int)((unsigned)src[0] << 8 | src[1]) << 16 | ((unsigned)src[2] << 8 | src[3]);
			src += 4;
		} while(--words);
	}

	// ==================
	// CSequenceUtilities
	// ==================

	const char* CSequenceUtilities::FOUR_BIT_PACKING   = "\xf\xf\xf\xf\xf\xf\xf\xf\xf\xf\xf\xf\xf\xf\xf\xf\xf\xf\xf\xf\xf\xf\xf\xf\xf\xf\xf\xf\xf\xf\xf\xf\xf\xf\xf\xf\xf\xf\xf\xf\xf\xf\xf\xf\xf\xc\xf\xf\xf\xf\xf\xf\xf\xf\xf\xf\xf\xf\xf\xf\xf\xf\xf\xf\xf\x0\xf\x1\xf\xf\xf\x2\xf\xf\xf\x9\xf\x4\xa\xf\xf\xf\x5\x7\x3\xf\xf\x6\xb\x8\xf";
	const char* CSequenceUtilities::FOUR_BIT_UNPACKING = "ACGTMRWSYKNX-XXX";

	// Java string hash function (32-bit)
	unsigned int CSequenceUtilities::JavaStringHash(const char* c, const unsigned int length) {

		if(length == 0) return 0;

		unsigned int hash = 0;

		const char* pSequence = c;
		for(unsigned int i = 0; i < length; i++) {
			hash = 31 * hash + *pSequence;
			pSequence++;
		}

		return hash;
	}

	// packs the reference and query bases into the specified packing string
	void CSequenceUtilities::Pack(string& packedString, const string& referenceBases, const string& queryBases) {

		// make sure that both strings are the same length
		if(referenceBases.size() != queryBases.size()) {
			printf("ERROR: Both strings must be the same length for 4-bit packing to succeed\n");
			exit(1);
		}

		// OBS: FOUR_BIT_PACKING translations have been checked

		const unsigned int packLength = queryBases.size();
		packedString.resize(packLength);
		char* pPackedString    = (char*)packedString.data();
		const char* pReference = (const char*)referenceBases.data();
		const char* pQuery     = (const char*)queryBases.data();

		// pack the current string
		for(unsigned int i = 0; i < packLength; i++) {
			pPackedString[i] = FOUR_BIT_PACKING[pReference[i]] | (FOUR_BIT_PACKING[pQuery[i]] << 4);
		}
	}

	// unpacks the reference and query bases from the specified packing string
	void CSequenceUtilities::Unpack(const string& packedString, string& referenceBases, string& queryBases) {

		const unsigned int packLength = packedString.size();
		referenceBases.resize(packLength);
		queryBases.resize(packLength);

		const char* pPackedString = (const char*)packedString.data();
		char* pReference          = (char*)referenceBases.data();
		char* pQuery              = (char*)queryBases.data();

		// OBS: FOUR_BIT_UNPACKING translations have been checked

		// unpack the current string
		for(unsigned int i = 0; i < packLength; i++) {		
			pQuery[i]     = FOUR_BIT_UNPACKING[(pPackedString[i] >> 4) & PACK_MASK];
			pReference[i] = FOUR_BIT_UNPACKING[pPackedString[i] & PACK_MASK];
		}
	}
}
