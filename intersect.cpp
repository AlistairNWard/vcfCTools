// ******************************************************
// vcfCTools (c) 2011 Alistair Ward
// Marth Lab, Department of Biology, Boston College
// All rights reserved.
// ------------------------------------------------------
// Last modified: 18 February 2011
// ------------------------------------------------------
// Calculate the intersection of two vcf files or a vcf
// file and a bed file.
// ******************************************************

#include "intersect.h"

using namespace std;
using namespace vcfCTools;

// Constructor.
intersect::intersect(void) {}

// Destructor.
intersect::~intersect(void) {}

// Set the Boolean flags required to determine the intersection operation to
// perform
void intersect::setBooleanFlags(bool findCommon, bool findUnion, bool findUnique, bool sitesOnly, bool annotate, bool within) {
  flags.findCommon   = findCommon;
  flags.findUnion    = findUnion;
  flags.findUnique   = findUnique;
  flags.sitesOnly    = sitesOnly;
  flags.annotate     = annotate;
  flags.whollyWithin = within;
}

// Intersect two vcf files.  Intersect by variant position only.
void intersect::intersectVcf(vcfHeader& header1, vcfHeader& header2, vcf& v1, variant& var1, vcf& v2, variant& var2, output& ofile) {
  bool write;

  // Build the variant structures for each vcf file.
  v1.success = v1.getRecord();
  v1.success = var1.buildVariantStructure(v1);
  v2.success = v2.getRecord();
  v2.success = var2.buildVariantStructure(v2);

  // Set the pointers to the start of each variant map.  For comparison
  // purposes, it is easiest to loop over variantMap structures.  These
  // structures contain all variants at a particular locus after any
  // excess bases have been trimmed from the reference/alternate alleles.
  // For output, the originalVariantsMap is used and needs to be cleared
  // to avoid using too much memory.  When all variants at a particular
  // locus in this map have been compared, the entry is sent to the output
  // and erased.
  var1.vmIter = var1.variantMap.begin();
  var2.vmIter = var2.variantMap.begin();
  var1.ovmIter = var1.originalVariantsMap.begin();
  var2.ovmIter = var2.originalVariantsMap.begin();

  // Parse and compare the two variant structures until the end of one of the files
  // is reached and the variant structure for that file is empty.
  while ( !(var1.variantMap.size() == 0 && !v1.success) || !(var2.variantMap.size() == 0 && !v2.success) ) {

    // If the two variant structures are built with the same reference sequence, compare
    // the contents and parse through all varians for this reference sequence.
    if (var1.vmIter->second.referenceSequence == var2.vmIter->second.referenceSequence) {
      string currentReferenceSequence = var1.vmIter->second.referenceSequence;

      // Since there are records from both vcf files containing this reference
      // sequence, set the reference sequence information variable, usedInComparison
      // to true.  If the vcf file is not sorted, the intersection will not work
      // correctly, but the vcf file should fail validation and, as such, shouldn't
      // be used.
      if (var1.variantMap.size() != 0) {var1.referenceSequenceInfo[currentReferenceSequence].usedInComparison = true;}
      if (var2.variantMap.size() != 0) {var2.referenceSequenceInfo[currentReferenceSequence].usedInComparison = true;}

      while (var1.variantMap.size() != 0 && var2.variantMap.size() != 0) {

        // Variants at the same locus.
        if (var1.vmIter->first == var2.vmIter->first) {
          var1.compareVariantsSameLocus(var2, flags);

          // Clear the compared variants from the structure and add the next one from 
          // the file into the structure if it is from the same reference sequence.
          var1.variantMap.erase(var1.vmIter);
          if (v1.variantRecord.referenceSequence == currentReferenceSequence && v1.success) {
            var1.addVariantToStructure(v1.position, v1.variantRecord);
            v1.success = v1.getRecord();
          }
          if (var1.variantMap.size() != 0) {var1.vmIter = var1.variantMap.begin();}

          var2.variantMap.erase(var2.vmIter);
          if (v2.variantRecord.referenceSequence == currentReferenceSequence && v2.success) {
            var2.addVariantToStructure(v2.position, v2.variantRecord);
            v2.success = v2.getRecord();
          }
          if (var2.variantMap.size() != 0) {var2.vmIter = var2.variantMap.begin();}

        // Variant from the first vcf file is at a larger coordinate than that in the
        // second vcf file.  Parse through the second file until the position is greater
        // than or equal to that in the second file.
        } else if (var1.vmIter->first > var2.vmIter->first) {
          if (flags.findCommon) {var2.filterUnique();}
          var2.variantMap.erase(var2.vmIter);
          if (v2.variantRecord.referenceSequence == currentReferenceSequence && v2.success) {
            var2.addVariantToStructure(v2.position, v2.variantRecord);
            v2.success = v2.getRecord();
          }

          // Reset the iterator if there are still variants in the structure.
          if (var2.variantMap.size() != 0) {
            var2.vmIter = var2.variantMap.begin();

          // If the variant position in the first file is greater than that in the
          // second and there are no more variants left for this reference sequence
          // in the second file, clear out the rest of the variants from the first
          // file for this reference sequence.
          } else {

            // First clear out any variants left in the originalVariants structure for
            // the second file.
            write = (flags.annotate) ? false : !flags.writeFromFirst;
            var2.clearOriginalVariants(header2, flags, ofile, write);

            // Then clear the remaining variants from the first file.
            if (var1.originalVariantsMap.size() != 0) {
              var1.clearReferenceSequence(header1, v1, flags, currentReferenceSequence, ofile, flags.writeFromFirst);
            }
          }
          
        // Variant from the first vcf file is at a smaller coordinate than that in the
        // second vcf file.
        } else if (var1.vmIter->first < var2.vmIter->first) {
          if (flags.findCommon && !flags.annotate) {var1.filterUnique();}
          var1.variantMap.erase(var1.vmIter);
          if (v1.variantRecord.referenceSequence == currentReferenceSequence && v1.success) {
            var1.addVariantToStructure(v1.position, v1.variantRecord);
            v1.success = v1.getRecord();
          }

          // Reset the iterator if there are still variants in the structure.
          if (var1.variantMap.size() != 0) {
            var1.vmIter = var1.variantMap.begin();

          // If the variant position in the second file is greater than that in the
          // first and there are no more variants left for this reference sequence
          // in the first file, clear out the rest of the variants from the second
          // file for this reference sequence.
          } else {

            // First clear out any variants left in the originalVariants structure for
            // the first file.
            var1.clearOriginalVariants(header1, flags, ofile, flags.writeFromFirst);

            // Then clear the remaining variants from the second file.
            if (var2.originalVariantsMap.size() != 0) {
              write = (flags.annotate) ? false : !flags.writeFromFirst;
              var2.clearReferenceSequence(header2, v2, flags, currentReferenceSequence, ofile, write);
            }
          }
        }

        // If the current position is beyond the max position in the originalVariantsMap
        // then all of the variants in this position have been compared and so it can
        // be sent to the output and erased.
        if (var1.variantMap.size() != 0 && var1.originalVariantsMap.size() != 0) {
          while (var1.vmIter->first > var1.ovmIter->second.begin()->maxPosition && var1.originalVariantsMap.size() != 0) {
            if (flags.writeFromFirst || flags.findUnion) {var1.buildOutputRecord(ofile, header1);}
            var1.originalVariantsMap.erase(var1.ovmIter);
            if (var1.originalVariantsMap.size() != 0) {var1.ovmIter = var1.originalVariantsMap.begin();}
          }
        }

        // If all of the records in the variantMap are exhausted there is nothing else to
        // compare.  If there are still records that haven't been sent to the output, these
        // should be sent now.
        if (var1.variantMap.size() == 0) {var1.clearOriginalVariants(header1, flags, ofile, flags.writeFromFirst);}

        // In order to avoid large memory usage, also clear out the originalVariantsMap
        // for var2.  If this is an annotation task, no records from the second vcf file
        // should be sent to the output.
        if (var2.variantMap.size() != 0 && var2.originalVariantsMap.size() != 0) {
          while (var2.vmIter->first > var2.ovmIter->second.begin()->maxPosition && var2.originalVariantsMap.size() != 0) {
            if (!flags.annotate && (!flags.writeFromFirst || flags.findUnion) ) {var2.buildOutputRecord(ofile, header2);}
            var2.originalVariantsMap.erase(var2.ovmIter);
            if (var2.originalVariantsMap.size() != 0) {var2.ovmIter = var2.originalVariantsMap.begin();}
          }
        }
        // If all of the records in the variantMap are exhausted there is nothing else to
        // compare.  If there are still records that haven't been sent to the output, these
        // should be sent now.
        write = (flags.annotate) ? false : !flags.writeFromFirst;
        if (var2.variantMap.size() == 0) {var2.clearOriginalVariants(header2, flags, ofile, write);}
      }

      // Having finished comparing, there may still be variants left from one of the two files.
      // Check that the two variant structures are empty and if not, finish processing the
      // remaining variants for this reference sequence.
      if (var1.originalVariantsMap.size() != 0) {
        var1.clearReferenceSequence(header1, v1, flags, currentReferenceSequence, ofile, flags.writeFromFirst);
      }
      if (var2.originalVariantsMap.size() != 0) {
        write = (flags.annotate) ? false : !flags.writeFromFirst;
        var2.clearReferenceSequence(header2, v2, flags, currentReferenceSequence, ofile, write);
      }

      // Now both variant maps are exhausted, so rebuild the maps with the variants from the
      // next reference sequence in the file.
      if (v1.success) {v1.success = var1.buildVariantStructure(v1);}
      if (v2.success) {v2.success = var2.buildVariantStructure(v2);}
 
      if (var1.variantMap.size() != 0) {var1.vmIter = var1.variantMap.begin();}
      if (var2.variantMap.size() != 0) {var2.vmIter = var2.variantMap.begin();}
      if (var1.originalVariantsMap.size() != 0) {var1.ovmIter = var1.originalVariantsMap.begin();}
      if (var2.originalVariantsMap.size() != 0) {var2.ovmIter = var2.originalVariantsMap.begin();}

    // If the variant structures are from different reference sequences, parse through the
    // second vcf file until the next reference sequence is found.  If finding the union
    // of the variants unique to the second vcf file, write them out.
    } else {
      if (var2.variantMap.size() != 0) {
        write = (flags.annotate) ? false : !flags.writeFromFirst;
        var2.clearReferenceSequence(header2, v2, flags, var2.vmIter->second.referenceSequence, ofile, write);
      }
      var2.buildVariantStructure(v2);
      if (var2.variantMap.size() != 0) {var2.vmIter = var2.variantMap.begin();}
    }
  }

  // If the variant structures are not empty, there was a problem and a warning is given.
  if (var1.variantMap.size() != 0 || (var2.variantMap.size() != 0 && !flags.annotate) ) {
    cerr << "WARNING: Not all records were flushed out of the variant structure." << endl;
  }

  // Flush the output buffer.
  ofile.flushOutputBuffer();
}

// Intersect a vcf file and a bed file.  It is assumed that the 
// two files are sorted by genomic coordinates and the reference
// sequences are in the same order.  Do not group together variants
// in common reference sequence.
void intersect::intersectVcfBed(vcfHeader& header, vcf& v, variant& var, bed& b, bedStructure& bs, output& ofile) {
  bs.lastBedIntervalEnd      = 0;
  unsigned int distanceToBed = 0;
  unsigned int leftDistance  = 0;
  unsigned int rightDistance = 0;

  // Build the variant structure for the vcf file.
  v.success   = v.getRecord();
  v.success   = var.buildVariantStructure(v);
  var.ovmIter = var.originalVariantsMap.begin();

  // Build the bed file structure and check if the first interval overlaps the
  // next.  If so, the bed intervals need to be broken up into unique
  // intervals.
  bs.lastBedInterval = false;
  b.success          = b.getRecord();
  bs.initialiseBedMap(b, flags);

  // Perform the intersection by looping over all records in the vcf file.  Keep
  // looping until this file is exhausted.  Within this loop, loop over each
  // reference sequence individually.
  while (v.success) {

    // Define the current reference sequence as that from the first entry in the
    // variant map.
    currentReferenceSequence = var.ovmIter->second.begin()->referenceSequence;

    // Loop over all records for this reference sequence.
    while (var.ovmIter->second.begin()->referenceSequence == currentReferenceSequence) {

      // Reset the iterateVcf and iterateBed flags.  Having compared all of the
      // variants appearing at this position, this will determine whether the
      // vcf file or the bed file should be iterated.
      iterateBed = false;
      iterateVcf = false;

      // Loop over all variants at this locus.
      var.ovIter = var.ovmIter->second.begin();
      for (; var.ovIter != var.ovmIter->second.end(); var.ovIter++) {

        // Consider each variant allele in turn and use the filtered vector to
        // indicate if a particular allele should be removed.  The actual removal
        // and modification of the genotypes (if they exist) are dealt with when
        // the variants are written to file.
        //
        // Begin by defining the required iterators.
        vector<int>::iterator posIter          = var.ovIter->reducedPosition.begin();
        vector<string>::iterator refIter       = var.ovIter->reducedRef.begin();
        vector<string>::iterator altIter       = var.ovIter->reducedAlts.begin();
        vector<variantType>::iterator typeIter = var.ovIter->type.begin();
        vector<bool>::iterator filtIter        = var.ovIter->filtered.begin();

        // Loop over the variants.
        for (; posIter != var.ovIter->reducedPosition.end(); posIter++) {

          // Determine the coordinate of the last base in the allele.  This
          // is used to determine if the alleles are wholly within the bed
          // interval.
          int endPos = max((refIter->size() + *posIter - 1), (altIter->size() + *posIter -1));

          // Define the start and end of the bed interval.
          int bedStart = bs.bmIter->first;
          int bedEnd   = bs.bmIter->second.end;

          // Define the variants overlap with the interval.
          bool beforeStart    = (endPos < bedStart ) ? true : false;
          bool overlapStart   = (*posIter < bedStart && endPos >= bedStart) ? true : false;
          bool within         = (*posIter >= bedStart && endPos <= bedEnd) ? true : false;
          bool overlapEnd     = (*posIter <= bedEnd && endPos > bedEnd) ? true : false;
          bool afterEnd       = (*posIter > bedEnd) ? true : false;

          // Variant is prior to the bed interval.  This depends on whether the
          // ref and alt alleles are required to fall wholly within the bed
          // interval.
          if (beforeStart || (overlapStart && flags.whollyWithin)) {
            *filtIter = priorToInterval(flags);

          // Variant is beyond the bed interval.
          } else if (afterEnd || (overlapEnd && flags.whollyWithin)) {
            beyondInterval();

          // Variant is within the bed interval
          } else if (within || ( (overlapStart || overlapEnd) && !flags.whollyWithin)) {
            *filtIter = withinInterval(flags);

          // If the else statement is reached there is an error in the code.
          } else {
            cerr << "ERROR: Problem with intersection algorithm detected." << endl;
            cerr << "The vcf record at " << var.ovIter->referenceSequence;
            cerr << ":" << var.ovmIter->first << " does not fall before, within or" << endl;
            cerr << "after the bed interval.  An algorithmic error must be present." << endl;
            cerr << "Program terminated." << endl;
            exit(1);
          }

          // Iterate the remaining iterators.
          refIter++;
          altIter++;
          typeIter++;
          filtIter++;
        }
      }

      // Iterate the either the vcf or the bed file.
      if (iterateVcf) {
        iterateVcfFile(header, v, var, ofile);

        // If the variant structure is empty, either the vcf file has been completely
        // parsed or all records for the particular reference sequence have been
        // exhausted.  If this is the case, build the new variant structure for the
        // new reference sequence and parse the bed file until this sequence is found.
        if (var.originalVariantsMap.size() == 0 && v.success) {nextReferenceSequence(v, var, b, bs);}

      // Iterate the bed file.
      } else if (iterateBed) {
        iterateBedFile(b, bs);
      } else {
        cerr << "ERROR: Neither the vcf or bed file are to be iterated." << endl;
        cerr << "This should not occur and indicates that an error exists" << endl;
        cerr << "in the vcf/bed intersection algorithm." << endl << endl;
        cerr << "Program terminated." << endl;
        exit(1);
      }

      // Check if the bed file is exhausted.  If the vcf file is completely parsed, the
      // program will naturally terminate, however, if the bed file has been exhausted,
      // the loops will continue.  Thus, if the bed map is empty and there are no more
      // intervals to be read in, finish parsing the vcf file without checking for
      // intersection with bed intervals.  If only records falling within bed intervals
      // are required, the program can terminate.
      if (bs.bedMap.size() == 0) {

        // The bed map is empty, but the bed file is not.  This means that all of the
        // intervals for the current reference sequence have been parsed.  The
        // remaining variants for this reference sequence can be processed since they
        // cannot overlap an interval and the next reference sequence can be started.
        if (b.success) {
          var.clearReferenceSequenceBed(header, v, flags, currentReferenceSequence, ofile);
          v.success   = var.buildVariantStructure(v);
          var.ovmIter = var.originalVariantsMap.begin();
          currentReferenceSequence = var.ovmIter->second.begin()->referenceSequence;
          bs.lastBedInterval = false;
          bs.initialiseBedMap(b, flags);

        // The bed file is completely parsed, so clear out the remaining variants in the
        // vcf file and terminate.
        } else {
          if (flags.findUnique || flags.annotate) {
            var.clearOriginalVariants(header, flags, ofile, true);
  
          // No more intersections will be founda and annotations aren't being performed, 
          // so no more records in the vcf file will be sent to the output so the loops 
          // can terminate.
          } else {
            v.success = false;
          }
        }
      }
    }
  }

  // If the variant structure is not empty, not all of the records were parsed.
  if (var.variantMap.size() != 0) {cerr << "WARNING: Not all records were flushed out of the variant structure." << endl;}

  // Flush the output buffer.
  ofile.flushOutputBuffer();
}

// When comparing a vcf with a bed, deal with the variants that fall before a bed
// interval.
bool intersect::priorToInterval(intFlags& flags) {
  // Update statistics on distance to the nearest bed interval.
  //unsigned int leftDistance  = *posIter - lastBedIntervalEnd;
  //unsigned int rightDistance = bedStart - *posIter;
  //if (lastBedIntervalEnd == 0) {distanceToBed = rightDistance;}
  //else {distanceToBed = (rightDistance > leftDistance) ? -1 * leftDistance : rightDistance;}
  //distanceDist[currentReferenceSequence][distanceToBed]++;

  // Mark this as a variant to be written out if variants unique to the
  // vcf file (i.e. outside of the bed intervals) were requested.
  bool filter  = flags.findUnique ? false : true;
  iterateVcf = true;

  return filter;
}

// When comparing a vcf with a bed, deal with the variants that fall after a bed
// interval.
void intersect::beyondInterval() {

  // Update statistics on distance to the nearest bed interval.
  //if (lastBedInterval) {
  //  distanceToBed = *posIter - bedEnd;
  //  distanceDist[currentReferenceSequence][distanceToBed]++;
  //}
  iterateBed = true;
}

// When comparing a vcf with a bed, deal with the variants that fall within a bed
// interval.
bool intersect::withinInterval(intFlags& flags) {

  // Update statistics on distance to the nearest bed interval.
  //distanceToBed = 0;
  //distanceDist[currentReferenceSequence][distanceToBed]++;

  // Mark this as a variant to be written out if variants intersectiong the
  // vcf file (i.e. within the bed intervals) were requested.
  bool filter  = flags.findUnique ? true : false;
  iterateVcf = true;

  return filter;
}

// After comparing all of the variants at a particular position, either
// the vcf or the bed file needs to be moved on one more record and
// relevant map entries can be erased.  Perform these tasks here.
//
// Begin with iterating the vcf file.
void intersect::iterateVcfFile(vcfHeader& header, vcf& v, variant& var, output& ofile) {

  // Build the output record, removing unwanted alleles and modifying the
  // genotypes if necessary and send to the output buffer.
  var.buildOutputRecord(ofile, header);

  // Erase the parsed and compared variant.
  var.originalVariantsMap.erase(var.ovmIter);

  // Add the next variant record from the current reference sequence into the structure.
  if (v.variantRecord.referenceSequence == currentReferenceSequence && v.success) {
    var.addVariantToStructure(v.position, v.variantRecord);
    v.success = v.getRecord();
  }

  // Set the iterator to the first element in the map.
  var.ovmIter = var.originalVariantsMap.begin();
}

// Or iterate the bed file.
void intersect::iterateBedFile(bed& b, bedStructure& bs) {
  bs.lastBedIntervalEnd = bs.bmIter->second.end;

  // Erase the last bed interval from the structure.
  bs.bedMap.erase(bs.bmIter);

  // Add the next bed interval from the current reference sequence into the structure.
  if (b.bRecord.referenceSequence == currentReferenceSequence && b.success) {
    bs.addIntervalToStructure(b.bRecord);
    b.success = b.getRecord();
  }

  // Set the iterator to the first element in the map.
  bs.bmIter = bs.bedMap.begin();

  // Set the bsNext iterator to the next element in the map and determine if this is
  // the last interval in the bed file for this reference sequence.
  bs.bmNext = bs.bedMap.begin();
  if (bs.bedMap.size() == 1) {bs.lastBedInterval = true;}
  else if (bs.bedMap.size() != 0) {bs.bmNext++;}

  // Check if the end position of the current interval is larger than the start
  // position of the next interval in the structure.  If so, the intervals
  // overlap and need to be split up.
  if (!bs.lastBedInterval) {
    if (bs.bmIter->second.end >= bs.bmNext->first) {
      bs.resolveOverlaps(flags.annotate);
      bs.bmIter = bs.bedMap.begin();
      bs.bmNext = bs.bedMap.begin();
      bs.bmNext++;
    }
  }
}

// Build the variant structure for the new reference sequence and move the bed
// file on until it reaches this reference sequence.
void intersect::nextReferenceSequence(vcf& v, variant& var, bed& b, bedStructure& bs) {
  v.success                = var.buildVariantStructure(v);
  var.ovmIter              = var.originalVariantsMap.begin();
  currentReferenceSequence = var.ovmIter->second.begin()->referenceSequence;

  // Clear the current bed map.
  bs.bmIter  = bs.bedMap.begin();
  for (; bs.bmIter != bs.bedMap.end(); bs.bmIter++) {bs.bedMap.erase(bs.bmIter);}

  // Parse through the bed file until the current reference sequence is found (or
  // the end of the bed file is reached).
  while (b.bRecord.referenceSequence != currentReferenceSequence) {b.success = b.getRecord();}
  if (b.success) {
    bs.lastBedInterval = false;
    bs.initialiseBedMap(b, flags);
  }
}

// Check to see that records for all reference sequences were compared
// correctly.
void intersect::checkReferenceSequences(variant& var1, variant& var2) {
  bool error = false;

  // Loop through all of the observed reference sequences from the first
  // variant structure and for each one, check if that reference sequence
  // exists in the second structure.  If so, check to see if the
  // usedInComparison is set to true.  If not, the ordering of the records
  // in one of the files is different to the other and so the variants
  // were not compared to each other.  This will lead to erroneous results.
  var1.refSeqIter = var1.referenceSequenceInfo.begin();
  var2.refSeqIter = var2.referenceSequenceInfo.begin();
  for (; var1.refSeqIter != var1.referenceSequenceInfo.end(); var1.refSeqIter++) {

    // Reference sequence not present in second variant structure.
    if (var2.referenceSequenceInfo.count(var1.refSeqIter->first) == 0) {

    // Reference sequence present in second variant structure.
    } else {
      if (!var2.refSeqIter->second.usedInComparison) {
        cerr << "WARNING: Variants were not compared for reference sequence: " << var2.refSeqIter->first << endl;
        error = true;
      }

      // If the records for this reference sequence were not contiguous, give a
      // warning.
      if (!var2.refSeqIter->second.contiguous) {
        cerr << "WARNING: Non-contiguous records in file b for reference sequence: " << var2.refSeqIter->first << endl;;
        error = true;
      }

      // Erase this element from the reference sequence information.
      var2.referenceSequenceInfo.erase(var2.refSeqIter);
      if (var2.referenceSequenceInfo.size() != 0) {var2.refSeqIter = var2.referenceSequenceInfo.begin();}
    }

    // If the records for this reference sequence were not contiguous, give a
    // warning.
    if (!var1.refSeqIter->second.contiguous) {
      cerr << "WARNING: Non-contiguous records in file a for reference sequence: " << var1.refSeqIter->first << endl;;
      error = true;
    }
  }

  // Display a final warning if errors were discovered.
  if (error) {
    cerr << endl;
    cerr << "WARNING: Not all reference sequences that should have been compared, were compared." << endl;
    cerr << "This is most likely due to an unsorted vcf file, or that the order in which the" << endl;
    cerr << "reference sequences appeared in the two files is different." << endl;
    cerr << endl;
    cerr << "Verify that the vcf files are sorted and reference sequences appear in the same" << endl;
    cerr << "order.  vcfCTools validate can be used to help establish this." << endl;
  }
}
