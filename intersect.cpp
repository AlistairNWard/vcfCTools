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
void intersect::intersectVcf(vcf& v1, variant& var1, vcf& v2, variant& var2, output& ofile) {

  // Build the variant structures for each vcf file.
  v1.success = var1.buildVariantStructure(v1);
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

  // The following Boolean flags are used when adding variants to the map.  They
  // determine whether or not to write out records that are unique to one of the
  // files.  write1 is set to true if the vcf file is being annotated.  This
  // ensures that every record from the vcf file is written out regardless of
  // whether it intersects with the other vcf file.  In the case that the two
  // files intersect, the additional annotation routine will be called.
  //bool write1 = ( (flags.findUnique && writeFrom == "a") || flags.findUnion || flags.annotate) ? true : false;
  //bool write2 = ( (flags.findUnique && writeFrom == "b") || flags.findUnion) ? true : false;

  // Parse and compare the two variant structures until the end of one of the files
  // is reached and the variant structure for that file is empty.
  while ( !(var1.variantMap.size() == 0 && !v1.success) && !(var2.variantMap.size() == 0 && !v2.success) ) {

    // If the two variant structures are built with the same reference sequence, compare
    // the contents and parse through all varians for this reference sequence.
    if (var1.vmIter->second.referenceSequence == var2.vmIter->second.referenceSequence) {
      string currentReferenceSequence = var1.vmIter->second.referenceSequence;

      while (var1.variantMap.size() != 0 && var2.variantMap.size() != 0) {

        // Variants at the same locus.
        if (var1.vmIter->first == var2.vmIter->first) {
          var1.compareVariantsSameLocus(var2, flags);

          // Clear the compared variants from the structure and add the next one from 
          // the file into the structure if it is from the same reference sequence.
          var1.variantMap.erase(var1.vmIter);
          if (v1.variantRecord.referenceSequence == currentReferenceSequence && v1.success) {
            var1.addVariantToStructure(v1.position, v1.variantRecord);
            v1.success = v1.getRecord(currentReferenceSequence);
          }
          if (var1.variantMap.size() != 0) {var1.vmIter = var1.variantMap.begin();}

          var2.variantMap.erase(var2.vmIter);
          if (v2.variantRecord.referenceSequence == currentReferenceSequence && v2.success) {
            var2.addVariantToStructure(v2.position, v2.variantRecord);
            v2.success = v2.getRecord(currentReferenceSequence);
          }
          if (var2.variantMap.size() != 0) {var2.vmIter = var2.variantMap.begin();}

        // Variant from the first vcf file is at a larger coordinate than that in the
        // second vcf file.  Parse through the second file until the position is greater
        // than or equal to that in the second file.
        } else if (var1.vmIter->first > var2.vmIter->first) {

//**************************************************
//**************************************************
//**************************************************
//**************************************************
// UNION NEEDS TO HAVE VARIANTS FROM THE SECOND FILE.
//**************************************************
//**************************************************
//**************************************************
//**************************************************

          var2.variantMap.erase(var2.vmIter);
          if (v2.variantRecord.referenceSequence == currentReferenceSequence && v2.success) {
            var2.addVariantToStructure(v2.position, v2.variantRecord);
            v2.success = v2.getRecord(currentReferenceSequence);
          }
          if (var2.variantMap.size() != 0) {var2.vmIter = var2.variantMap.begin();}
          
        // Variant from the first vcf file is at a smaller coordinate than that in the
        // second vcf file.
        } else if (var1.vmIter->first < var2.vmIter->first) {
          if (!flags.findUnique && !flags.annotate) {var1.filterUnique();}
          var1.variantMap.erase(var1.vmIter);
          if (v1.variantRecord.referenceSequence == currentReferenceSequence && v1.success) {
            var1.addVariantToStructure(v1.position, v1.variantRecord);
            v1.success = v1.getRecord(currentReferenceSequence);
          }
          if (var1.variantMap.size() != 0) {var1.vmIter = var1.variantMap.begin();}
        }

        // If the current position is beyond the max position in the originalVariantsMap
        // then all of the variants in this position have been compared and so it can
        // be sent to the output and erased.
        while (var1.vmIter->first > var1.ovmIter->second.maxPosition && var1.originalVariantsMap.size() != 0) {
          var1.buildOutputRecord(ofile);
          var1.originalVariantsMap.erase(var1.ovmIter);
          if (var1.originalVariantsMap.size() != 0) {var1.ovmIter = var1.originalVariantsMap.begin();}
        }
      }

      // Having finished comparing, there may still be variants left from one of the two files.
      // Check that the two variant structures are empty and if not, finish processing the
      // remaining variants for this reference sequence.
      if (var1.originalVariantsMap.size() != 0) {var1.clearReferenceSequence(v1, flags, currentReferenceSequence, ofile, true);}
      if (var2.originalVariantsMap.size() != 0) {var2.clearReferenceSequence(v2, flags, currentReferenceSequence, ofile, false);}

      // Now both variant maps are exhausted, so rebuild the maps with the variants from the
      // next reference sequence in the file.
      v1.success = var1.buildVariantStructure(v1);
      v2.success = var2.buildVariantStructure(v2);
 
      if (var1.variantMap.size() != 0) {var1.vmIter = var1.variantMap.begin();}
      if (var2.variantMap.size() != 0) {var2.vmIter = var2.variantMap.begin();}

    // If the variant structures are from different reference sequences, parse through the
    // second vcf file until the next reference sequence is found.  If finding the union
    // of the variants unique to the second vcf file, write them out.
    } else {
      if (var2.variantMap.size() != 0) {var2.clearReferenceSequence(v2, flags, var2.vmIter->second.referenceSequence, ofile, false);}
      var2.buildVariantStructure(v2);
      if (var2.variantMap.size() != 0) {var2.vmIter = var2.variantMap.begin();}
    }
  }

  // If the variant structures are not empty, there was a problem and a warning is given.
  if (var1.variantMap.size() != 0 || var2.variantMap.size() != 0) {
    cerr << "WARNING: Not all records were flushed out of the variant structure." << endl;
  }

  // Flush the output buffer.
  ofile.flushOutputBuffer();
}

// Intersect a vcf file and a bed file.  It is assumed that the 
// two files are sorted by genomic coordinates and the reference
// sequences are in the same order.  Do not group together variants
// in common reference sequence.
void intersect::intersectVcfBed(vcf& v, variant& var, bed& b, bedStructure& bs, output& ofile) {
  map<int, bedRecord>::iterator bmNext;
  unsigned int lastBedIntervalEnd = 0;
  unsigned int distanceToBed      = 0;
  unsigned int leftDistance       = 0;
  unsigned int rightDistance      = 0;
  bool lastBedInterval            = false;

  // Build the variant structures for the vcf and bed files.
  v.success = var.buildVariantStructure(v);
  b.success = bs.buildBedStructure(b);

  // Set the pointers to the start of each variant map.  For intersection with
  // a bed file, there is no need to work with the reduced allele structures,
  // so the original variants structure is used.
  var.ovmIter = var.originalVariantsMap.begin();
  bs.bmIter   = bs.bedMap.begin();
  bmNext      = bs.bedMap.begin();
  if (bs.bedMap.size() == 1) {
    lastBedInterval = true;
  } else {
    bmNext++;
  }

  // Check if the end position of the current interval is larger than the start
  // position of the next interval in the structure.  If so, the intervals
  // overlap and need to be split up.
  if (!lastBedInterval) {
    if (bs.bmIter->second.end >= bmNext->first) {
      bs.resolveOverlaps(flags.annotate);
      bs.bmIter  = bs.bedMap.begin();
      bmNext     = bs.bedMap.begin();
      bmNext++;
    }
  }

  // Define the current reference sequence as that from the first entry in the
  // variant map.
  string currentReferenceSequence = var.ovmIter->second.referenceSequence;

// Parse through the vcf file until it is finished and the variant structure is empty.
  //while (var.variantMap.size() != 0 && bs.bedMap.size() != 0) {
  while (var.originalVariantsMap.size() != 0 || v.success) {
    if (var.ovmIter->second.referenceSequence == bs.bmIter->second.referenceSequence) {

      // Consider each variant allele in turn and use the filtered vector to
      // indicate if a particular allele should be removed.  The actual removal
      // and modification of the genotypes (if they exist) are dealt with when
      // the variants are written to file.
      //
      // Begin by defining the required iterators.
      vector<int>::iterator posIter          = var.ovmIter->second.reducedPosition.begin();
      vector<string>::iterator refIter       = var.ovmIter->second.reducedRef.begin();
      vector<string>::iterator altIter       = var.ovmIter->second.reducedAlts.begin();
      vector<variantType>::iterator typeIter = var.ovmIter->second.type.begin();
      vector<bool>::iterator filtIter        = var.ovmIter->second.filtered.begin();

      // Initialise flags that indicate whether to iterate the bed file
      // or the vcf file.
      bool iterateBed = false;
      bool iterateVcf = false;

      // Loop over the variants.
      for (; posIter != var.ovmIter->second.reducedPosition.end(); posIter++) {

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
          leftDistance  = *posIter - lastBedIntervalEnd;
          rightDistance = bedStart - *posIter;
          //distanceToBed = (lastBedIntervalEnd == 0) ? rightDistance : min(leftDistance, rightDistance);
          if (lastBedIntervalEnd == 0) {
            distanceToBed = rightDistance;
          } else {
            distanceToBed = (rightDistance > leftDistance) ? -1 * leftDistance : rightDistance;
          }
          distanceDist[currentReferenceSequence][distanceToBed]++;

          // Mark this as a variant to be written out if variants unique to the
          // vcf file (i.e. outside of the bed intervals) were requested.
          *filtIter  = flags.findUnique ? false : true;
          iterateVcf = true;

        // Variant is beyond the bed interval.
        } else if (afterEnd || (overlapEnd && flags.whollyWithin)) {
          if (lastBedInterval) {
            distanceToBed = *posIter - bedEnd;
            distanceDist[currentReferenceSequence][distanceToBed]++;
          }
          iterateBed = true;

        // Variant is within the bed interval
        } else if (within || ( (overlapStart || overlapEnd) && !flags.whollyWithin)) {
          distanceToBed = 0;
          distanceDist[currentReferenceSequence][distanceToBed]++;
          *filtIter  = flags.findUnique ? true : false;
          iterateVcf = true;

        // If the else statement is reached there is an error in the code.
        } else {
          cerr << "ERROR: Problem with intersection algorithm detected." << endl;
          cerr << "The vcf record at " << var.ovmIter->second.referenceSequence;
          cerr << ":" << var.ovmIter->first << " does not fall before, within or" << endl;
          cerr << "after the bed interval.  An algorithmic error must be present." << endl;
          cerr << "Program terminated." << endl;
          exit(1);
        }

        // If annotations are required, perform them here.
        if (flags.annotate) {}

        // Iterate the remaining iterators.
        refIter++;
        altIter++;
        typeIter++;
        filtIter++;
      }

      // Iterate whichever file is necessary.  For example, if all variants were
      // prior to the bed interval, the vcf file will be interated until an
      // overlap with the bed file is found.
      if (iterateVcf) {

        // Build the output record, removing unwanted alleles and modifying the
        // genotypes if necessary and send to the output buffer.
        var.buildOutputRecord(ofile);

        // Update the originalVariants structure.
        var.originalVariantsMap.erase(var.ovmIter);
        if (v.variantRecord.referenceSequence == currentReferenceSequence && v.success) {
          var.addVariantToStructure(v.position, v.variantRecord);
          v.success = v.getRecord(currentReferenceSequence);
        }
        if (var.originalVariantsMap.size() != 0) {var.ovmIter = var.originalVariantsMap.begin();}

      // Iterate the bed file.
      } else if (iterateBed) {
        lastBedIntervalEnd = bs.bmIter->second.end;
        bs.bedMap.erase(bs.bmIter);
        if (b.bRecord.referenceSequence == currentReferenceSequence && b.success) {
          bs.addIntervalToStructure(b.bRecord);
          b.success = b.getRecord();
        }
        if (bs.bedMap.size() != 0) {
          bs.bmIter = bs.bedMap.begin();

          // This is the last bed interval if the size of the bed map is 1.
          bmNext = bs.bedMap.begin();
          if (bs.bedMap.size() == 1) {
            lastBedInterval = true;
          } else {
            bmNext++;
          }

          // Check if the end position of the current interval is larger than the start
          // position of the next interval in the structure.  If so, the intervals
          // overlap and need to be split up.
          if (!lastBedInterval) {
            if (bs.bmIter->second.end >= bmNext->first) {
              bs.resolveOverlaps(flags.annotate);
              bs.bmIter = bs.bedMap.begin();
              bmNext    = bs.bedMap.begin();
              bmNext++;
            }
          }

        // If the bed map is exhausted, clear out the variant map and move to the next
        // reference sequence.
        } else {

          // If there are variants in the variant map, clear them out and set the current
          // reference sequence to the next reference sequence in the vcf file.
          if (var.originalVariantsMap.size() != 0) {
            var.clearReferenceSequenceBed(v, flags, var.ovmIter->second.referenceSequence, ofile);
          }
          v.success = var.buildVariantStructure(v);
          if (var.originalVariantsMap.size() != 0) {var.ovmIter = var.originalVariantsMap.begin();}

          lastBedInterval = false;

          // Build the new bed structure.
          if (b.success) {
            b.success = bs.buildBedStructure(b);
            if (bs.bedMap.size() != 0) {bs.bmIter = bs.bedMap.begin();}
          }
        }

      // If neither requires iterating then something strange has happened and the
      // code is not behaving as expected.
      } else {
        cerr << "ERROR: Neither the vcf or bed file are to be iterated." << endl;
        cerr << "This should not occur and indicates that an error exists" << endl;
        cerr << "in the vcf/bed intersection algorithm." << endl << endl;
        cerr << "Program terminated." << endl;
        exit(1);
      }
    }
  }

// If the variant structure is not empty, not all of the records were parsed.
  if (var.variantMap.size() != 0) {cerr << "WARNING: Not all records were flushed out of the variant structure." << endl;}

  // Flush the output buffer.
  ofile.flushOutputBuffer();
}
