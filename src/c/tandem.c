#include <stdlib.h>

#include <stdio.h>   // for testing
#include <string.h>  // for testing

static const int MAX_DISPLACEMENT = 12;

unsigned char* twobit( char* sequence, int length, int offset ) {
  
  // size in bytes, rounded up to whole number of 64-bit words; and allowing for the maximum displacement
  int buflen = (1 + (((length + MAX_DISPLACEMENT) * 2 ) | 63)) / 8;
  unsigned char* buffer = malloc( (size_t)buflen );
  int current_nuc_index = offset;
  int i,j;
  for (i=0; i<buflen; i++) { // loop over bytes
    unsigned char byte = 0;
    for (j=0; j<8; j+=2) { // loop over bit position within byte
      switch (sequence[current_nuc_index]) {
      case 'A':
	break;
      case 'C':
	byte |= (1 << j);
	break;
      case 'G':
	byte |= (2 << j);
	break;
      case 'T':
	byte |= (3 << j);
	break;
      case 0:
	--current_nuc_index; // pad with As
	break;
      }
      ++current_nuc_index;
    }
    buffer[i] = byte;
  }
  return buffer;
}

inline void foundmatch( char* sizes, char* displacements, int pos, int size, int displacement, int length ) {
  if (pos + displacement + size > length) {
    size = length - displacement - pos;
  }
  // convert size to length of repetitive region
  size += displacement;
  // only accept true tandems
  if (size < 2*displacement)
    return;
  if (sizes[pos] < size) {
    sizes[pos] = size;
    displacements[pos] = displacement;
  }
}

inline int min(int i, int j) {
  if (i<j) return i;
  return j;
}

void annotate( char* sequence, char* sizes, char* displacements, int length) {

  unsigned char* seqs[4] = { twobit( sequence, length, 0 ), 
			     twobit( sequence, length, 1 ),
			     twobit( sequence, length, 2 ),
			     twobit( sequence, length, 3 ) };

  // initialize size and displacement arrays
  int pos, displacement;
  for (pos=0; pos < length; pos++) {
    sizes[pos] = 1;
    displacements[pos] = 1;
  }
  
  // loop over starting positions within the sequence
  for (pos=0; pos < length; pos+=4) {

    // get 128 bits, 64 nucleotides, at pos, in two 64-bit longs
    unsigned long long original0 = *((long long*) & seqs[0][pos/4]);
    unsigned long long original1 = *((long long*) & seqs[0][(pos+32)/4]);

    for (displacement=1; displacement < MAX_DISPLACEMENT; displacement+=1) {

      if (pos + displacement >= length) break;
      
      // get target
      unsigned long long target0 = *((long long*) & seqs[displacement % 4][(pos+displacement)/4]);

      // calculate match lengths
      target0 ^= original0;
      int size0 = 64;
      int size1 = 64;
      int size2 = 64;
      int size3 = 64;
      int nucscanleft = (__builtin_ffsll( target0 ) + 1)/2;  // 0 for no mismatches; 1 for nuc 0 mismatch, etc.
      if (nucscanleft == 1) {
	size0 = 0;                                         // found result for position 0
	target0 &= -4LL;                                   // remove blocking mismatch
	nucscanleft = (__builtin_ffsll( target0 ) + 1)/2;  // recompute
      }
      if (nucscanleft == 2) {
	size0 = min(size0,1);
	size1 = 0;
	target0 &= -16LL;
	nucscanleft = (__builtin_ffsll( target0 ) + 1)/2;
      }
      if (nucscanleft == 3) {
	size0 = min(size0, 2);
	size1 = min(size1, 1);
	size2 = 0;
	target0 &= -64LL;
	nucscanleft = (__builtin_ffsll( target0 ) + 1)/2;
      }
      if (nucscanleft == 0) {
	unsigned long long target1 = *((long long*) & seqs[displacement % 4][(pos+displacement+32)/4]);
	target1 ^= original1;
	nucscanleft = (__builtin_ffsll( target1 ) +1)/2;
	if (nucscanleft == 0) {
	  nucscanleft = 32;
	} else {
	  nucscanleft--;
	}
	size0 = min(size0, nucscanleft+32);
	size1 = min(size1, nucscanleft+31);
	size2 = min(size2, nucscanleft+30);
	size3 = nucscanleft + 29;
      } else {
	size0 = min(size0, nucscanleft-1);
	size1 = min(size1, nucscanleft-2);
	size2 = min(size2, nucscanleft-3);
	size3 = nucscanleft - 4;
      }
      foundmatch( sizes, displacements, pos, size0, displacement, length );
      foundmatch( sizes, displacements, pos+1, size1, displacement, length );
      foundmatch( sizes, displacements, pos+2, size2, displacement, length );
      foundmatch( sizes, displacements, pos+3, size3, displacement, length );
    }
  }
  free(seqs[0]);
  free(seqs[1]);
  free(seqs[2]);
  free(seqs[3]);
}

int main() {

  char* seq1 = "TATTTGCATGCGCTTTCGAGCTGTTGAAGAGACGTGTATTGGAATAAGTAATCACATAAGTGTTAGTAACTTATTTAAATACGTATAGAGTCGCCTATTTGCCTAGCCTTTTGGTTCTCAGATTTTTTAATTATTACATTGCTATAAGGGTGTAACTGTGTGATAGCCAAAATTTTAAGCTGCAAATGGTTTGTAAATATGATATATTACAAGCTTCATGAAAATCGGTTTATGACTGATCCGCGATTACGTTGAAAGGCGACTGGCAGAGATACTTTTGTTCAGATGTTTTTTCAGGTAGCGATTCCAATGAATAGGTAAAATACCTTGCAAGTTTTGTTGTTGTCGTTGGAGGAAATGTGGATGTGGTTGTTATTGTTGA";
  char* sizes = (char*)malloc( strlen(seq1)+1 );
  char* displacements = (char*)malloc( strlen(seq1)+1 );
  
  annotate( seq1, sizes, displacements, strlen(seq1) );

  int i;
  for (i=0; i<strlen(seq1); i++) {
    printf ("%c\t%u\t%u\n",seq1[i], sizes[i], displacements[i]);
  }
  return 0;

}
