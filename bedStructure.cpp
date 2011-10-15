// ******************************************************
// vcfCTools (c) 2011 Alistair Ward
// Marth Lab, Department of Biology, Boston College
// All rights reserved.
// ------------------------------------------------------
// Last modified: 18 February 2011
// ------------------------------------------------------
// Define the bedStructure class.  This is essentially the
// equivalent of the variant class, but for bed files and
// as such is significantly simpler.
// ******************************************************

#include "bedStructure.h"

using namespace std;
using namespace vcfCTools;

//Constructor.
bedStructure::bedStructure(void) {
  recordsInMemory = 1000;
};

// Destructor.
bedStructure::~bedStructure(void) {};

// Build up a new map containing the bed intervals.  Set the iterator to
// the first interval and check if it overlaps with the following interval.
void bedStructure::initialiseBedMap(bed& b, intFlags& flags) {
  buildBedStructure(b);

  // Set the pointers to the start of each variant map.  For intersection with
  // a bed file, there is no need to work with the reduced allele structures,
  // so the original variants structure is used.
  bmIter = bedMap.begin();
  bmNext = bedMap.begin();
  if (bedMap.size() == 1) {lastBedInterval = true;}
  else {bmNext++;}

  // Check if the end position of the current interval is larger than the start
  // position of the next interval in the structure.  If so, the intervals
  // overlap and need to be split up.
  if (!lastBedInterval) {
    if (bmIter->second.end >= bmNext->first) {
      resolveOverlaps(flags.annotate);
      bmIter = bedMap.begin();
      bmNext = bedMap.begin();
      bmNext++;
    }
  }
}

// Build up a structure containing bed intervals and annotations.
bool bedStructure::buildBedStructure(bed& b) {
  unsigned int count = 0;

  string tempReferenceSequence = b.bRecord.referenceSequence;
  while (b.success && count < recordsInMemory && b.bRecord.referenceSequence == tempReferenceSequence) {
    addIntervalToStructure(b.bRecord);
    b.success = b.getRecord();
    count++;
  }

  return b.success;
}

// Add the bed record to the structure.
void bedStructure::addIntervalToStructure(bedRecord& br) {

// If this position already has an entry, determine how to join the intervals.
  if (bedMap.count(br.start) > 0) {
    resolveOverlaps(br);
  } else {
    bedMap[br.start].referenceSequence = br.referenceSequence;
    bedMap[br.start].start = br.start;
    bedMap[br.start].end = br.end;
    bedMap[br.start].info = br.info;
  }
}

// If two bed records overlap, split them into two (or three if one of the
// records is wholly within another) records.  One of these records
// is a common interval and the info strings for each record need to be
// unpacked, compared for common strings and repacked with unique strings.
// The other record(s) is(are) the non-overlapping section and will retain
// the info string of its parent.
//
// The first of these routines is called when a record to be added has the
// same start coordinate of an existing record.
void bedStructure::resolveOverlaps(bedRecord& b1) {
  bedRecord b2 = bedMap[b1.start];

  // Erase the map element with the duplicated start position.
  bedMap.erase(b1.start);

  // The first interval is shorter than the second.
  if (b1.end < b2.end) {

    // The interval start - b1.end is common.  b1.end + 1 - b2.end is unique to b2.
    // Unique interval.
    bedMap[b1.end + 1].referenceSequence = b2.referenceSequence;
    bedMap[b1.end + 1].start = b1.end + 1;
    bedMap[b1.end + 1].end = b2.end;
    bedMap[b1.end + 1].info = b2.info;

    // Common interval.
    generateCommonInterval(b1, b2, b1.start, b1.end);

  // The two intervals are identical.
  } else if (b1.end == b2.end) {
    generateCommonInterval(b1, b2, b1.start, b1.end);

  // The second interval is shorter than the first.
  } else if (b1.end > b2.end) {

    // The interval start - b2.end is common.  b2.end + 1 - b1.end is unique to b1.
    // Unique interval.
    bedMap[b2.end + 1].referenceSequence = b1.referenceSequence;
    bedMap[b2.end + 1].start = b2.end + 1;
    bedMap[b2.end + 1].end = b1.end;
    bedMap[b2.end + 1].info = b1.info;

    // Common interval.
    generateCommonInterval(b1, b2, b1.start, b2.end);
  }
}

// Resolve overlap between two records already in the structure (e.g. they
// do not have the same start coordinate).
void bedStructure::resolveOverlaps(bool annotate) {
  bool resolvedIntervals = false;
  map<int, bedRecord>::iterator bmNext = bedMap.begin();
  map<int, int> startPositions, endPositions;
  map<int, int>::iterator sIter, eIter;
  vector<bedRecord> bRecords;
  vector<bedRecord>::iterator bIter;
  bmNext++;

  // Store the first bed record in the vector.
  bRecords.push_back(bmIter->second);
  startPositions[bmIter->first] = 1;
  endPositions[bmIter->second.end] = 1;

  // Determine how many consecutive records overlap.
  while (!resolvedIntervals && bmNext != bedMap.end()) {
    if (bmIter->second.end >= bmNext->first) {

      // Add the start and end coordinates to the maps.
      startPositions[bmNext->first] = 1;
      endPositions[bmNext->second.end] = 1;

      // Store the record in the vector, then delete from the map.  This
      // reocrd is no longer required in the map as it will be broken up.
      bRecords.push_back(bmNext->second);
      bedMap.erase(bmNext);

      // Get the second element of the map to check if this also overlaps
      // with the first record.
      bmNext = bedMap.begin();
      bmNext++;
    } else {
      resolvedIntervals = true;
    }
  }

  // Also remove the first element from the map.  This is also to be broken
  // up.
  bedMap.erase(bmIter);

  // Now build up the new records.  Proceed by selecting the first (leftmost
  // coordinate) start position and take this as the start of the first record.
  // The end coordinate is the smaller of the leftmost end coordinate (of any
  // of the stored records) or the start coordinate - 1.  Remove these elements
  // from the vectors.
  vector<int> reconstructedStarts, reconstructedEnds;

  // Start the reconstructed intervals with the start coordinate of the first
  // original bed record.
  sIter = startPositions.begin();
  reconstructedStarts.push_back(sIter->first);
  startPositions.erase(sIter);
  while (startPositions.size() != 0) {

    // Find the end position.
    sIter = startPositions.begin();
    eIter = endPositions.begin();

    // lowest (start - 1 ) < lowest( end )
    if (sIter->first - 1 < eIter->first) {
      reconstructedEnds.push_back(sIter->first - 1);

      // The start of the next reconstructed interval is this end coordinate + 1.
      reconstructedStarts.push_back(sIter->first);
      startPositions.erase(sIter);
      sIter = startPositions.begin();

    // lowest (start - 1 ) > lowest( end )
    } else {
      reconstructedEnds.push_back(eIter->first);
      
      // The start of the next reconstructed interval is this end coordinate + 1.
      reconstructedStarts.push_back(eIter->first + 1);
      endPositions.erase(eIter);
      eIter = endPositions.begin();
    }
  }

  // There will still be end positions to process.  Finish up the intervals
  // using these end positions.
  int end;
  while (endPositions.size() != 0) {
    eIter = endPositions.begin();
    end = eIter->first;
    reconstructedEnds.push_back(end);
    endPositions.erase(eIter);
    if (endPositions.size() != 0) {reconstructedStarts.push_back(end + 1);}
  }

  // Add these bed intervals back into the bed map.  Annotations will be dealt
  // with afterwards (if necessary).
  if (reconstructedStarts.size() != reconstructedEnds.size()) {
    cerr << "ERROR:  Problem in reconstructing overlapping bed records." << endl;
    cerr << "First record in region at: " << bRecords[0].referenceSequence << ":" << reconstructedStarts[0] << endl;
    exit(1);
  }
  vector<int>::iterator stIter, enIter;
  enIter = reconstructedEnds.begin();
  for (stIter = reconstructedStarts.begin(); stIter != reconstructedStarts.end(); stIter++) {
    bedMap[*stIter].referenceSequence = bRecords[0].referenceSequence;
    bedMap[*stIter].start = *stIter;
    bedMap[*stIter].end = *enIter;
    enIter++;
  }

  // Now that the coordinates of the reconstructed bed intervals are known,
  // determine the annotations that fall within these regions.
  if (annotate) {
    vector<string> ann;
    vector<string>::iterator annIter;
    map<string, int> annMap;
    map<string, int>::iterator annMapIter;
    string constructInfo;

    enIter = reconstructedEnds.begin();
    for (stIter = reconstructedStarts.begin(); stIter != reconstructedStarts.end(); stIter++) {
      constructInfo = "";
      bIter = bRecords.begin();
      while (bIter != bRecords.end()) {

        // If the end coordinate of the original record is less than the start
        // coordinate of this reconstructed record, the original record can
        // be deleted as it no longer contributes to the reconstructed
        // intervals.
        if (bIter->end < *stIter) {
          bIter = bRecords.erase(bIter);

        // This record overlaps the reconstructed interval, so add the annotation
        // information to the map.
        } else if (bIter->start <= *enIter) {
          ann = split(bIter->info, ";");
          for (annIter = ann.begin(); annIter != ann.end(); annIter++) {
            annMap[*annIter] = 1;
          }
          bIter++;
        } else {
          bIter++;
        }
      }

      // Add the info string to the bed map and clear the annotation map.
      for (annMapIter = annMap.begin(); annMapIter != annMap.end(); annMapIter++) {
        constructInfo = (constructInfo == "") ? annMapIter->first : constructInfo + ";" + annMapIter->first;
        annMap.erase(annMapIter);
      }
      bedMap[*stIter].info = constructInfo;

      enIter++;
    }
  }
}

// Generate the entry for the common interval.
void bedStructure::generateCommonInterval(bedRecord& b1, bedRecord& b2, int start, int end) {
  string info = "";
  map<string, int> existingInfo;

  // Unpack the info strings from the two intersecting records.
  vector<string> b1Info = split(b1.info, ";");
  vector<string> b2Info = split(b2.info, ";");

  // Generate the new interval.
  bedMap[start].referenceSequence = b1.referenceSequence;
  bedMap[start].start = start;
  bedMap[start].end = end;
  for (vector<string>::iterator i = b1Info.begin(); i != b1Info.end(); i++) {
    info = (info == "") ? *i : info + ";" + *i;
    existingInfo[*i] = 0;
  }
  for (vector<string>::iterator i = b2Info.begin(); i != b2Info.end(); i++) {
    if (existingInfo.count(*i) == 0) {info = (info == "") ? *i : info + ";" + *i;}
  }
  bedMap[start].info = info;
}
