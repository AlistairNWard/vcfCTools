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

// Intersect two vcf files.  Intersect by variant position only.
void intersectVcf(vcf& v1, variant& var1, vcf& v2, variant& var2, bool findUnion, bool findUnique, bool annotate, string writeFrom, ostream* output) {

  // The following Boolean flags are used when adding variants to the map.  They
  // determine whether or not to write out records that are unique to one of the
  // files.  write1 is set to true if the vcf file is being annotated.  This
  // ensures that every record from the vcf file is written out regardless of
  // whether it intersects with the other vcf file.  In the case that the two
  // files intersect, the additional annotation routine will be called.
  bool write1 = ( (findUnique && writeFrom == "a") || findUnion || annotate) ? true : false;
  bool write2 = ( (findUnique && writeFrom == "b") || findUnion) ? true : false;

  // Build the variant structures for each vcf file.
  v1.success = var1.buildVariantStructure(v1);
  v2.success = var2.buildVariantStructure(v2);

  // Set the pointers to the start of each variant map.
  var1.vmIter = var1.variantMap.begin();
  var2.vmIter = var2.variantMap.begin();

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
          if (!findUnique) {
            if (annotate) {var1.annotateRecordVcf(var2.vmIter->second, v2.dbsnpVcf);}
            var1.writeVariants(output);
          }

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
          if (write2) {var2.writeVariants(output);}
          var2.variantMap.erase(var2.vmIter);
          if (v2.variantRecord.referenceSequence == currentReferenceSequence && v2.success) {
            var2.addVariantToStructure(v2.position, v2.variantRecord);
            v2.success = v2.getRecord(currentReferenceSequence);
          }
          if (var2.variantMap.size() != 0) {var2.vmIter = var2.variantMap.begin();}
          
        // Variant from the first vcf file is at a smaller coordinate than that in the
        // second vcf file.
        } else if (var1.vmIter->first < var2.vmIter->first) {
          if (write1) {var1.writeVariants(output);}
          var1.variantMap.erase(var1.vmIter);
          if (v1.variantRecord.referenceSequence == currentReferenceSequence && v1.success) {
            var1.addVariantToStructure(v1.position, v1.variantRecord);
            v1.success = v1.getRecord(currentReferenceSequence);
          }
          if (var1.variantMap.size() != 0) {var1.vmIter = var1.variantMap.begin();}
        }
      }

      // Having finished comparing, there may still be variants left from one of the two files.
      // Check that the two variant structures are empty and if not, finish processing the
      // remaining variants for this reference sequence.
      if (var1.variantMap.size() != 0) {var1.clearReferenceSequence(v1, currentReferenceSequence, write1, output);}
      if (var2.variantMap.size() != 0) {var2.clearReferenceSequence(v2, currentReferenceSequence, write2, output);}

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
      if (var2.variantMap.size() != 0) {var2.clearReferenceSequence(v2, var2.vmIter->second.referenceSequence, write2, output);}
      var2.buildVariantStructure(v2);
      if (var2.variantMap.size() != 0) {var2.vmIter = var2.variantMap.begin();}
    }
  }
}

// Intersect a vcf file and a bed file.  It is assumed that the 
// two files are sorted by genomic coordinates and the reference
// sequences are in the same order.  Do not group together variants
// in common reference sequence.
void intersectVcfBed(vcf& v, variant& var, bed& b, bedStructure& bs, bool findUnique, bool annotate, ostream* output) {
  map<int, bedRecord>::iterator bmNext;
  bool lastBedInterval = false;

  // Build the variant structures for the vcf and bed files.
  v.success = var.buildVariantStructure(v);
  b.success = bs.buildBedStructure(b);

  // Set the pointers to the start of each variant map.
  var.vmIter = var.variantMap.begin();
  bs.bmIter = bs.bedMap.begin();
  bmNext = bs.bedMap.begin();
  bmNext++;
  if (bmNext == bs.bedMap.end()) {lastBedInterval = true;}

  // Check if the end position of the current interval is larger than the start
  // position of the next interval in the structure.  If so, the intervals
  // overlap and need to be split up.
  if (!lastBedInterval) {
    if (bs.bmIter->second.end >= bmNext->first) {bs.resolveOverlaps(annotate);}
  }

  string currentReferenceSequence = v.referenceSequence;

// As soon as the end of the vcf file is reached, there are no
// more intersections and the program can terminate.
  while (v.success && var.variantMap.size() != 0 && bs.bedMap.size() != 0) {
    if (var.vmIter->second.referenceSequence == bs.bmIter->second.referenceSequence) {

      // Variant is prior to the bed interval.
      if (var.vmIter->first < bs.bmIter->first) {
        if (findUnique || annotate) {var.writeVariants(output);}
        var.variantMap.erase(var.vmIter);
        if (v.variantRecord.referenceSequence == currentReferenceSequence && v.success) {
          var.addVariantToStructure(v.position, v.variantRecord);
          v.success = v.getRecord(currentReferenceSequence);
        }
        if (var.variantMap.size() != 0) {var.vmIter = var.variantMap.begin();}

      // Variant has start coordinate past the bed interval.
      } else if (var.vmIter->first > bs.bmIter->second.end) {
        bs.bedMap.erase(bs.bmIter);
        if (b.bRecord.referenceSequence == currentReferenceSequence && b.success) {
          bs.addIntervalToStructure(b.bRecord);
          b.success = b.getRecord();
        }
        if (bs.bedMap.size() != 0) {
          bs.bmIter = bs.bedMap.begin();
          if (bs.bedMap.begin()++ != bs.bedMap.end()) {
            bmNext = bs.bedMap.begin();
            bmNext++;
          } else {lastBedInterval = true;}

          // Check if the end position of the current interval is larger than the start
          // position of the next interval in the structure.  If so, the intervals
          // overlap and need to be split up.
          if (!lastBedInterval) {
            if (bs.bmIter->second.end >= bmNext->first) {bs.resolveOverlaps(annotate);}
          }

        // If the bed map is exhausted, clear out the variant map and move to the next
        // reference sequence.
        } else {

          // If there are variants in the variant map, clear them out and set the current
          // reference sequence to the next reference sequence in the vcf file.
          if (var.variantMap.size() != 0) {
            bool write = (findUnique || annotate) ? true : false;
            var.clearReferenceSequence(v, var.vmIter->second.referenceSequence, write, output);
          }
          v.success = var.buildVariantStructure(v);
          if (var.variantMap.size() != 0) {var.vmIter = var.variantMap.begin();}

          lastBedInterval = false;

          // Build the new bed structure.
          b.success = bs.buildBedStructure(b);
          if (bs.bedMap.size() != 0) {bs.bmIter = bs.bedMap.begin();}
        }

      // Variant lies within the bed interval.
      } else {
        if (annotate) {var.annotateRecordBed(bs.bmIter->second);}
        if (!findUnique) {var.writeVariants(output);}
        var.variantMap.erase(var.vmIter);
        if (v.variantRecord.referenceSequence == currentReferenceSequence && v.success) {
          var.addVariantToStructure(v.position, v.variantRecord);
          v.success = v.getRecord(currentReferenceSequence);
        }
        if (var.variantMap.size() != 0) {var.vmIter = var.variantMap.begin();}
      }

    // If the reference sequence of the bed interval is not the same as that of the vcf
    // record,
    } else {

      // If there are variants in the variant map, clear them out and set the current
      // reference sequence to the next reference sequence in the vcf file.
      if (var.variantMap.size() != 0) {
        bool write = (findUnique || annotate) ? false : true;
        var.clearReferenceSequence(v, var.vmIter->second.referenceSequence, write, output);
      }
      v.success = var.buildVariantStructure(v);
      if (var.variantMap.size() != 0) {var.vmIter = var.variantMap.begin();}

      lastBedInterval = false;

      // Clear the bed map.
      bs.bedMap.clear();

      // Read through the remaining records until the correct reference sequence is
      // found.  There is no need to put the records in the map, since the vcf will
      // not be intersected with these records.
      while (b.bRecord.referenceSequence != currentReferenceSequence && b.success) {b.success = b.getRecord();}

      // Build the new bed structure.
      b.success = bs.buildBedStructure(b);
      if (bs.bedMap.size() != 0) {bs.bmIter = bs.bedMap.begin();}
    }
  }
}
