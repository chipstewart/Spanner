// ***************************************************************************
// MosaikAlignment.h - a concatenation of MOSAIK classes that allows the user
//                     to read from and write to MOSAIK alignment archives.
// ---------------------------------------------------------------------------
// (c) 2006 - 2009 Michael Strömberg
// Marth Lab, Department of Biology, Boston College
// ---------------------------------------------------------------------------
// Dual licenced under the GNU General Public License 2.0+ license or as
// a commercial license with the Marth Lab.
// ***************************************************************************

#pragma once

#include <cmath>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <iostream>
#include <map>
#include <string>
#include <vector>
#include "fastlz.h"

#ifndef WIN32
#include <tr1/unordered_map>
#else
#include <unordered_map>
#include "windows.h"
#endif

#ifndef __APPLE__
#include <cstdint>
#endif

using namespace std;
using namespace std::tr1;

#ifdef WIN32
typedef unsigned long long uint64_t;
#endif

namespace Mosaik {

	// ==========
	// data types
	// ==========

#define SIZEOF_OFF_TYPE    8

#ifdef WIN32
#define ftell64(a)     _ftelli64(a)
#define fseek64(a,b,c) _fseeki64(a,b,c)
	typedef __int64 off_type;
#elif SPARC
#define ftell64(a)     ftello(a)
#define fseek64(a,b,c) fseeko(a,b,c)
	typedef off64_t off_type;
#elif __APPLE__
#define ftell64(a)     ftello(a)
#define fseek64(a,b,c) fseeko(a,b,c)
#define fstat64(a,b)   fstat(a,b)
#define stat64 stat
	typedef off_t off_type;
#else
#define ftell64(a)     ftello(a)
#define fseek64(a,b,c) fseeko(a,b,c)
	typedef off_t off_type;
#endif

#define SIZEOF_CHAR          1
#define SIZEOF_WCHAR         2
#define SIZEOF_SHORT         2
#define SIZEOF_INT           4
#define SIZEOF_FLOAT         4
#define SIZEOF_DOUBLE        8
#define SIZEOF_UINT64        8

#define FASTLZ_BETTER_SPEED       1
#define FASTLZ_BETTER_COMPRESSION 2

#define NUM_READS_OFFSET 25

	// =========
	// Alignment
	// =========

	struct Alignment {
		unsigned int MateReferenceBegin;   // required for SAM/BAM
		unsigned int MateReferenceEnd;     // required for SAM/BAM
		unsigned int MateReferenceIndex;   // required for SAM/BAM
		unsigned int ReferenceBegin;
		unsigned int ReferenceEnd;
		unsigned int ReferenceIndex;
		unsigned int ReadGroupCode;        // the read group code (temp)
		unsigned short NumMismatches;      // number of mismatches
		unsigned short QueryBegin;
		unsigned short QueryEnd;
		unsigned char Quality;             // alignment quality
		bool IsFirstMate;                  // is this alignment from the first mate in a paired-end read
		bool IsMateReverseStrand;          // read orientation for the mate
		bool IsPairedEnd;                  // is the read sequenced as a paired-end read
		bool IsResolvedAsPair;             // is the alignment part of resolved paired-end read
		bool IsReverseStrand;              // read orientation
		bool WasRescued;                   // was the alignment rescued during local alignment search
		char* ReferenceName;               // only filled via CAlignmentReader (temp)
		string Reference;
		string Query;
		string BaseQualities;
		string Name;                       // the read name

		// constructors
		Alignment(void)
			: MateReferenceBegin(0)
			, MateReferenceEnd(0)
			, MateReferenceIndex(0)
			, ReferenceIndex(0)
			, ReadGroupCode(0)
			, Quality(0)
			, IsFirstMate(false)
			, IsMateReverseStrand(false)
			, IsPairedEnd(false)
			, IsResolvedAsPair(false)
			, IsReverseStrand(false)
			, WasRescued(false)
			, ReferenceName(NULL)
		{}

		// our less-than operator
		bool operator<(const Alignment& al) const {
			if(ReferenceIndex == al.ReferenceIndex) return ReferenceBegin < al.ReferenceBegin;
			return ReferenceIndex < al.ReferenceIndex;
		}
	};

	// ===============
	// AlignmentStatus
	// ===============

	// define our alignment status flags (relevant before MosaikMerge)
	typedef unsigned char AlignmentStatus;
#define AS_UNKNOWN                      0   // specifies an unset status flag
#define AS_SINGLE_END_READ              1   // transferred from the read format
#define AS_PAIRED_END_READ              2   // transferred from the read format
#define AS_UNSORTED_READ                4   // expected in MosaikAligner data
#define AS_SORTED_ALIGNMENT             8   // expected in MosaikSort data
#define AS_ALL_MODE                     16  // enables non-unique PE resolution
#define AS_UNIQUE_MODE                  32  // disables non-unique PE resolution

	// define our read flags
#define RF_UNKNOWN                      0   // specifies an unset read flag
#define RF_FAILED_QUALITY_CHECK         1   // reserved, not currently used
#define RF_HAVE_MATE1                   2   // specifies if mate1 is stored
#define RF_HAVE_MATE2                   4   // specifies if mate2 is stored
#define RF_IS_LONG_READ                 8   // specifies any read that is longer than 255 bases
#define RF_IS_PAIRED_IN_SEQUENCING      16  // originates from a paired-end data set
#define RF_IS_PCR_OR_OPTICAL_DUPLICATE  32  // reserved, not currently used
#define RF_IS_UNALIGNED                 64  // reserved, not currently used
#define RF_RESOLVED_AS_PAIR             128 // seen only after MosaikSort or MosaikMerge

	// define our alignment flags
#define AF_UNKNOWN                      0   // specifies an unset alignment flag
#define AF_IS_FIRST_MATE                1   // reserved, not currently used
#define AF_IS_MATE_REVERSE_STRAND       2   // specifies the orientation of the mate
#define AF_IS_NOT_PRIMARY               4   // reserved, not currently used
#define AF_IS_REVERSE_STRAND            8   // specifies the orientation of the current alignment
#define AF_IS_SECOND_MATE               16  // reserved, not currently used
#define AF_WAS_RESCUED                  32  // specifies if the alignment was rescued during local alignment search

	// define our alignment tags
#define AT_UNKNOWN                      0

	// define our header tags
#define HT_UNKNOWN                      0

	// define our reference sequence tags
#define RST_UNKNOWN                     0

	// define our read group tags
#define RGT_UNKNOWN                     0

	// define our tag data structure
	typedef unsigned char TagType;

	const TagType TT_CHAR   = 'c';
	const TagType TT_DOUBLE = 'd';
	const TagType TT_FLOAT  = 'f';
	const TagType TT_INT16  = 's';
	const TagType TT_INT32  = 'i';
	const TagType TT_INT64  = 'l';
	const TagType TT_STRING = 'z';
	const TagType TT_UCHAR  = 'C';
	const TagType TT_UINT16 = 'S';
	const TagType TT_UINT32 = 'I';
	const TagType TT_UINT64 = 'L';

#define TAG_STRING_LEN 512

	struct Tag {
		unsigned char ID;
		TagType Type;
		union {
			char Char;
			char String[TAG_STRING_LEN];
			double Double;
			float Float;
			int Int32;
			long long Int64;
			short Int16;
			unsigned char UChar;
			unsigned int UInt32;
			unsigned long long UInt64;
			unsigned short UInt16;
		};
	};

	// ===========
	// AlignedRead
	// ===========

	struct AlignedRead {
		unsigned int ReadGroupCode;
		string Name;
		vector<Alignment> Mate1Alignments;
		vector<Alignment> Mate2Alignments;
		bool IsLongRead;
		bool IsPairedEnd;

		// constructor
		AlignedRead()
			: ReadGroupCode(0)
			, IsLongRead(false)
			, IsPairedEnd(false)
		{}
	};

	// ======================
	// SequencingTechnologies
	// ======================

	typedef unsigned short SequencingTechnologies;

#define ST_UNKNOWN               0
#define ST_454                   1
#define ST_HELICOS               2
#define ST_ILLUMINA              4
#define ST_PACIFIC_BIOSCIENCES   8
#define ST_SOLID                16
#define ST_SANGER               32

	// =====
	// CSHA1
	// =====

#define NUM_SHA_IV_WORDS     5
#define PGP_SHA_BLOCKBYTES  64
#define NUM_SHA_KEY_WORDS   16
#define SHA_HASH_LENGTH     20
#define SHA_KEY_LENGTH      32
#define MAX_STRING_LENGTH  512

	// The SHA f()-functions
#define f1(x,y,z) (z ^ (x & (y ^ z)))			// Rounds  0 - 19
#define f2(x,y,z) (x ^ y ^ z)					// Rounds 20 - 39
#define f3(x,y,z) ((x & y) + (z & (x ^ y) ))	// Rounds 40 - 59
#define f4(x,y,z) (x ^ y ^ z)					// Rounds 60 - 79

	// The SHA Mysterious Constants.
#define K2	0x5A827999L	// Rounds  0-19 - floor(sqrt(2) * 2^30)
#define K3	0x6ED9EBA1L	// Rounds 20-39 - floor(sqrt(3) * 2^30)
#define K5	0x8F1BBCDCL	// Rounds 40-59 - floor(sqrt(5) * 2^30)
#define K10	0xCA62C1D6L	// Rounds 60-79 - floor(sqrt(10) * 2^30)

	// 32-bit rotate left
#define ROTL(n,X) ((X << n) | (X >> (32-n)))

	// The prototype SHA sub-round
#define subRound(a, b, c, d, e, f, k, data) \
	(e += ROTL(5,a) + f(b, c, d) + k + data, b = ROTL(30, b))

	// The initial expanding function
#define expandx(W,i) (t = W[i&15] ^ W[(i-14)&15] ^ W[(i-8)&15] ^ W[(i-3)&15], ROTL(1, t))
#define expand(W,i) (W[i&15] = expandx(W,i))

#define swap_int32(x) \
	(((x & 0x000000ff) << 24) + \
	((x & 0x0000ff00) <<  8) + \
	((x & 0x00ff0000) >>  8) + \
	((x & 0xff000000) >> 24))

	class CSHA1 {
	public:
		// generates the read group code based on the supplied text
		static unsigned int GenerateReadGroupCode(const string& readGroupID, const string& sampleName);

	private:
		// define our SHA context
		struct SHA1_Data {
			unsigned int key[NUM_SHA_KEY_WORDS];
			unsigned int iv[NUM_SHA_IV_WORDS];
			unsigned int bytes;
		};
		// final wrapup - pad to 64-byte boundary
		static void FinalizeSHA1(SHA1_Data& data, unsigned int* key, const unsigned char numWords);
		// initialize the SHA values
		static void InitializeSHA1(SHA1_Data& data);
		// Shuffle the bytes into big-endian order within words, as per the SHA spec.
		static void SwapBytes(unsigned int* dest, unsigned char const *src, unsigned char words);
		// perform the SHA transformation
		static void TransformSHA1(unsigned int* block, unsigned int* key);
		// update SHA for a block of data
		static void UpdateSHA1(SHA1_Data& data, void const *bufIn, size_t len);
		// our passphrase key
		unsigned char mKey[SHA_KEY_LENGTH];
	};

	// ==================
	// CSequenceUtilities
	// ==================

#define PACK_MASK     15

	class CSequenceUtilities {
	public:
		// Java string hash function (32-bit)
		static unsigned int JavaStringHash(const char* c, const unsigned int length);
		// packs the reference and query bases into the specified packing string
		static void Pack(string& packedString, const string& referenceBases, const string& queryBases);
		// unpacks the reference and query bases from the specified packing string
		static void Unpack(const string& packedString, string& referenceBases, string& queryBases);
	private:
		// the packing vector
		static const char* FOUR_BIT_PACKING;
		// the unpacking vector
		static const char* FOUR_BIT_UNPACKING;
	};

	// =========
	// ReadGroup
	// =========

	struct ReadGroup {
		unsigned int MedianFragmentLength;
		unsigned int ReadGroupCode;
		SequencingTechnologies SequencingTechnology;
		string CenterName;
		string Description;
		string LibraryName;
		string PlatformUnit;
		string ReadGroupID;
		string SampleName;

		// constructor
		ReadGroup(void)
			: MedianFragmentLength(0)
			, ReadGroupCode(0)
			, SequencingTechnology(ST_UNKNOWN)
		{}

		// create a 32-bit identifier for the read group code
		static unsigned int GetCode(const ReadGroup& readGroup) {
			const unsigned int readGroupCode = CSHA1::GenerateReadGroupCode(readGroup.ReadGroupID, readGroup.SampleName);
			return readGroupCode;return readGroupCode;
		}
	};

	typedef vector<ReadGroup> ReadGroupVector;

	// =================
	// ReferenceSequence
	// =================

	struct ReferenceSequence {
		uint64_t NumAligned;
		unsigned int Begin;
		unsigned int End;
		unsigned int NumBases;
		string Name;
		string Bases;
		string GenomeAssemblyID;
		string Species;
		string MD5;
		string URI;

		// constructor
		ReferenceSequence()
			: NumAligned(0)
			, Begin(0)
			, End(0)
			, NumBases(0)
		{}
	};

	typedef vector<ReferenceSequence> RefVector;

	// ================
	// CAlignmentReader
	// ================

	class CAlignmentReader {
	public:
		// constructor
		CAlignmentReader(void);
		// destructor
		~CAlignmentReader(void);
		// checks to see if this is truly a MOSAIK alignment archive
		static bool CheckFile(const string& filename, const bool showError);
		// closes the alignment archive
		void Close(void);
		// returns the number of bases in the archive
		uint64_t GetNumBases(void) const;
		// returns the number of reads in the archive
		uint64_t GetNumReads(void) const;
		// retrieves the read group data given a read group code
		ReadGroup GetReadGroupFromCode(const unsigned int code);
		// retrieves the read groups vector
		ReadGroupVector GetReadGroups(void) const;
		// retrieves the reference sequence data
		RefVector GetReferenceData(void) const;
		// retrieves the file status
		AlignmentStatus GetStatus(void) const;
		// jumps to the block containing the specified reference index and position
		void Jump(const unsigned int referenceIndex, const unsigned int referencePosition);
		// loads the next read from the alignment archive
		bool LoadNextRead(Mosaik::AlignedRead& ar);
		// opens the alignment archive
		void Open(const string& filename);
		// sets the file pointer to the beginning of the read data
		void Rewind(void);

	private:
		// load the read header from disk
		void LoadReadHeader(string& readName, unsigned int& readGroupCode, unsigned char& readStatus, unsigned int& numMate1Alignments, unsigned int& numMate2Alignments);
		// deserializes each alignment and stores them in the supplied vector
		void ReadAlignments(vector<Alignment>& alignments, const bool isLongRead, const bool isPairedInSequencing, const bool isResolvedAsPair, const unsigned int readGroupCode);
		// deserialize the alignment
		void ReadAlignment(Alignment& al, const bool isLongRead, const bool isPairedInSequencing, const bool isResolvedAsPair);
		// reads a new compressed partition (returns false if EOF occurs)
		bool ReadPartition(void);
		// reads the tag from disk
		void ReadTag(Tag& tag);
		// denotes the status of the output stream
		bool mIsOpen;
		// our compressed output stream
		FILE* mInStream;
		// stores the archive read count
		uint64_t mNumReads;
		uint64_t mNumBases;
		// stores the current read number
		uint64_t mCurrentRead;
		// stores the file offsets
		off_type mReadsOffset;
		off_type mReferenceGapOffset;
		off_type mIndexOffset;
		// our input buffer
		char* mBuffer;
		char* mBufferPtr;
		unsigned int mBufferLen;
		// our input compression buffer
		unsigned char* mCompressionBuffer;
		unsigned int mCompressionBufferLen;
		// our input filename
		string mInputFilename;
		// our partitioning setup
		unsigned short mPartitionSize;
		unsigned short mPartitionMembers;
		// our reference sequence LUT
		char** mRefSeqLUT;
		unsigned short mNumRefSeqs;
		// our reference sequences
		RefVector mReferenceSequences;
		// our read groups
		vector<ReadGroup> mReadGroups;
		// our file status
		AlignmentStatus mStatus;
		// our sequencing technology
		SequencingTechnologies mSeqTech;
		// our header tags
		map<unsigned char, Tag> mHeaderTags;
		// our read group LUT
		map<unsigned int, ReadGroup> mReadGroupLUT;
	};

	// ================
	// CAlignmentWriter
	// ================

#define MEMORY_BUFFER_SIZE 3000

	class CAlignmentWriter {
	public:
		// constructor
		CAlignmentWriter(void);
		// destructor
		~CAlignmentWriter(void);
		// adds a header tag
		void AddHeaderTag(const Tag& tag);
		// closes the alignment archive
		void Close(void);
		// retrieves the number of bases written
		uint64_t GetNumBases(void) const;
		// retrieves the number of reads written
		uint64_t GetNumReads(void) const;
		// opens the alignment archive
		void Open(const string& filename, const RefVector& referenceSequences, const vector<ReadGroup>& readGroups, const AlignmentStatus as);
		// saves the read to the alignment archive
		void SaveAlignedRead(const Mosaik::AlignedRead& ar);
		// adds a header tag (only works before opening the file)
		void AddHeaderTag(const unsigned char tagID, const TagType& tagType);

	private:
		// specifies our index entry
		struct IndexEntry {
			off_type Offset;
			unsigned int Position;
			unsigned int ReferenceIndex;
		};
		// adjusts the buffer
		void AdjustBuffer(void);
		// serializes the specified alignment
		void WriteAlignment(const Alignment* pAl, const bool isLongRead, const bool isPairedEnd, const bool isFirstMate, const bool isResolvedAsPair);
		// write partition to disk
		void WritePartition(void);
		// write the read header to disk
		void WriteReadHeader(const string& readName, const unsigned int readGroupCode, const unsigned char readStatus, const unsigned int numMate1Alignments, const unsigned int numMate2Alignments);
		// writes the tag to disk
		void WriteTag(const map<unsigned char, Tag>::const_iterator& htIter);
		// denotes the status of the output stream
		bool mIsOpen;
		// our compressed output stream
		FILE* mOutStream;
		// stores the number of sequences that have been written
		uint64_t mNumReads;
		uint64_t mNumBases;
		// our output buffer
		unsigned char* mBuffer;
		unsigned int mBufferLen;
		unsigned int mBufferPosition;
		unsigned int mBufferThreshold;
		// our output compression buffer
		unsigned char* mCompressionBuffer;
		unsigned int mCompressionBufferLen;
		// our partitioning setup
		unsigned short mPartitionSize;
		unsigned short mPartitionMembers;
		// our output filename
		string mOutputFilename;
		// our reference sequence vector
		RefVector mReferenceSequences;
		unsigned int mNumRefSeqs;
		// our reference sequence gap vector
		vector<unordered_map<unsigned int, unsigned short> >* mpRefGapVector;
		// our alignment status
		AlignmentStatus mStatus;
		// denotes that this alignment archive is paired-end (used in SaveRead)
		bool mIsPairedEndArchive;
		// our block index
		vector<IndexEntry> mIndex;
		unsigned int mLastReferenceIndex;
		unsigned int mLastReferencePosition;
		bool mStoreIndex;
		// our header tags
		map<unsigned char, Tag> mHeaderTags;
	};

	// ================
	// CMemoryUtilities
	// ================

	class CMemoryUtilities {
	public:
		// checks if the buffer is large enough to accommodate the requested size
		static void CheckBufferSize(char* &pBuffer, unsigned int& bufferLen, const unsigned int requestedBytes);
		// checks if the buffer is large enough to accommodate the requested size
		static void CheckBufferSize(unsigned char* &pBuffer, unsigned int& bufferLen, const unsigned int requestedBytes);
	};

	// ============
	// CTimeSupport
	// ============

	class CTimeSupport {
	public:
		// returns the current time
		static uint64_t GetSystemTime(void);
	private:
		// converts a time_t variable to our 64-bit notation
		static uint64_t ConvertTimeT(const time_t& time);
		// constants
		static const uint64_t SECS_BETWEEN_EPOCHS;
		static const uint64_t SECS_TO_100NS;
	};

	// =========
	// CFastLZIO
	// =========

	// the buffer is currently set to 10 MB
#define FASTLZ_IO_BUFFER_LEN        10485760

	// the buffer must be at least 5 % larger
#define FASTLZ_IO_OUTPUT_BUFFER_LEN  9986438

	class CFastLZIO {
	public:
		// constructor
		CFastLZIO(void);
		// destructor
		~CFastLZIO(void);
		// clears the buffer
		void Clear(void);
		// our input method
		void Read(char* &buffer, unsigned int& bufferLen, FILE* stm);
		// our input method (STL string)
		void Read(string& s, FILE* stm);
		// our output method
		void Write(const char* buffer, const unsigned int bufferLen, FILE* stm);
	private:
		// our buffer
		char* mBuffer;
		unsigned int mBufferLen;
	};
}

