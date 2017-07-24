#include "align.h"

/*****************************************************************************************************************
 This code is copyright (c) Gerton Lunter, Dec 2009
 It may not be distributed, made public, or used in other software without the permission of the copyright holder
******************************************************************************************************************/

#include <emmintrin.h>
#include <stdio.h>

// These are defined for use in the alignment routines
#define match_label 0
#define insert_label 1
#define delete_label 3

//_________________________________________________________________________________________________

int fastAlignmentRoutine(char* seq1, char* seq2, char* qual2, int len1, int len2, int gapextend, int nucprior, char* homopolgapq_or_localgapopen, char usehomopolgapq)
{
    // setup stuff
    // seq2 is the read; the shorter of the sequences

    // the bottom-left and top-right corners of the DP table are just
    // included at the extreme ends of the diagonal, which measures
    // n=8 entries diagonally across.  This fixes the length of the
    // longer (horizontal, seq1) sequence to 14 (2*8-2) more than the shorter

    const short pos_inf = 0x2000;
    const short gap_extend = gapextend*4;
    const short nuc_prior = nucprior*4;
    const short n_score = 1*4;

    short _score = 0;
    short minscore = pos_inf;
    short minscoreidx = -1;

    register __m128i _m1;
    register __m128i _i1;
    register __m128i _d1;
    register __m128i _m2;
    register __m128i _i2;
    register __m128i _d2;

    __m128i _seq1win;
    __m128i _seq2win;
    __m128i _qual2win;
    __m128i _seq1nqual;   // 0 if N, 0xffff if not
    __m128i _gap_extend = _mm_set1_epi16( gap_extend );
    __m128i _nuc_prior = _mm_set1_epi16( nuc_prior );
    __m128i _initmask = _mm_set_epi16( 0,0,0,0,0,0,0,-1 );

    // initialization
    _m1 = _mm_set_epi16( pos_inf, pos_inf, pos_inf, pos_inf, pos_inf, pos_inf, pos_inf, pos_inf );
    _i1 = _m1;
    _d1 = _m1;
    _m2 = _m1;
    _i2 = _m1;
    _d2 = _m1;
    _seq1nqual = _m1;

    // sequence 1 is initialized with the n-long prefix, in forward direction
    // sequence 2 is initialized as empty; reverse direction
    _seq1win = _mm_set_epi16( seq1[7], seq1[6], seq1[5], seq1[4], seq1[3], seq1[2], seq1[1], seq1[0] );
    _seq2win = _m1;
    _qual2win = _mm_set1_epi16(64*4);

    // if N, make n_score; if != N, make n_score+pos_inf
    _seq1nqual = _mm_and_si128( _mm_add_epi16( _mm_cmpeq_epi16( _seq1win, _mm_set1_epi16( 'N' ) ), _mm_set1_epi16( pos_inf + n_score ) ), _mm_set1_epi16( pos_inf + n_score ) );

    short inith[len1];
    int idx = len1;

    if (usehomopolgapq) {

      int homopol = -1;
      int homopollen = 0;

      while (idx) {
	  --idx;
	  if (seq1[idx] == homopol) {
	    homopollen += !!homopolgapq_or_localgapopen[homopollen+1];
	  } else {
	    homopollen = 0;
	  }
	  inith[idx] = 4*(homopolgapq_or_localgapopen[homopollen] - '!');
	  homopol = seq1[idx];
	  if (homopol == 'N') {
	    homopol=0;
	  }
      }
    } else {
      while (idx) {
	--idx;
	inith[idx] = homopolgapq_or_localgapopen[idx]*4;
      }
    }

    //debug -- print local gap open penalties
    //for (idx=0; idx<len1; idx++) printf("%c",((inith[idx]/4)+33) ); printf("\n");
    __m128i _gap_open = _mm_set_epi16(inith[7],inith[6],inith[5],inith[4],inith[3],inith[2],inith[1],inith[0]);

    // main loop.  Do one extra iteration, with nucs from sequence 2 just moved out
    // of the seq2win/qual arrays, to simplify getting back pointers
    int s = 0;

    for (; s < 2*(len2+8-1 + 1); s += 2) {

        const int sOver2 = s/2;

        // seq1 is current; seq2 needs updating
        _seq2win = _mm_slli_si128( _seq2win, 2 );
        _qual2win = _mm_slli_si128( _qual2win, 2 );

        if (sOver2 < len2) {
            _seq2win = _mm_insert_epi16( _seq2win, seq2[sOver2], 0 );
            _qual2win = _mm_insert_epi16( _qual2win, 4*(qual2[sOver2]), 0 );

        } else {
            _seq2win = _mm_insert_epi16( _seq2win, '0', 0 );
            _qual2win = _mm_insert_epi16( _qual2win, 64*4, 0 );
        }

        //
        // S even
        //

        _m1 = _mm_andnot_si128( _initmask, _m1 );
        _m2 = _mm_andnot_si128( _initmask, _m2 );
        _m1 = _mm_min_epi16( _m1, _mm_min_epi16( _i1, _d1 ) );

        // at this point, extract minimum score.  Referred-to position must
        // be y==len2-1, so that current position has y==len2; i==0 so d=0 and y=s/2

        if (sOver2 >= len2)
        {
            switch (sOver2 - len2)
            {
                case 0: _score = _mm_extract_epi16( _m1, 0 ); break;
                case 1: _score = _mm_extract_epi16( _m1, 1 ); break;
                case 2: _score = _mm_extract_epi16( _m1, 2 ); break;
                case 3: _score = _mm_extract_epi16( _m1, 3 ); break;
                case 4: _score = _mm_extract_epi16( _m1, 4 ); break;
                case 5: _score = _mm_extract_epi16( _m1, 5 ); break;
                case 6: _score = _mm_extract_epi16( _m1, 6 ); break;
                case 7: _score = _mm_extract_epi16( _m1, 7 ); break;

                default:
                    printf("Something is wrong in alignment switch statement\n");
            }

            if (_score < minscore)
            {
                minscore = _score;
                minscoreidx = s;     // point back to the match state at this entry, so as not to
            }                      // have to store the state at s-2
        }

        _m1 = _mm_add_epi16( _m1,
                _mm_min_epi16( _mm_andnot_si128( _mm_cmpeq_epi16( _seq2win,
                            _seq1win ),
                        _qual2win ),
                    _seq1nqual ) );
        _d1 = _mm_min_epi16( _mm_add_epi16( _d2,
                    _gap_extend ),
                _mm_add_epi16( _mm_min_epi16( _m2,
                        _i2 ),   // allow I->D
                    _gap_open ) );
        _d1 = _mm_insert_epi16( _mm_slli_si128( _d1, 2 ), pos_inf, 0 );
        _i1 = _mm_add_epi16( _mm_min_epi16( _mm_add_epi16( _i2,
                        _gap_extend ),
                    _mm_add_epi16( _m2,
                        _gap_open ) ),
                _nuc_prior );

        //
        // S odd
        //

        // seq1 needs updating; seq2 is current
        const char c = (8 + sOver2 < len1) ? seq1[ 8+(sOver2) ] : 'N';

        _seq1win = _mm_insert_epi16( _mm_srli_si128( _seq1win, 2 ), c, 8-1 );
        _seq1nqual = _mm_insert_epi16( _mm_srli_si128( _seq1nqual, 2 ), (c=='N')?n_score:pos_inf, 8-1 );

        _gap_open = _mm_insert_epi16( _mm_slli_si128( _gap_open, 2 ),
                inith[ 8 + sOver2 < len1 ? 8 + sOver2 : len1-1 ],
                0 );

        _initmask = _mm_slli_si128( _initmask, 2 );
        _m2 = _mm_min_epi16( _m2, _mm_min_epi16( _i2, _d2 ) );

        // at this point, extract minimum score.  Referred-to position must
        // be y==len2-1, so that current position has y==len2; i==0 so d=0 and y=s/2

        if (sOver2 >= len2)
        {
            switch (sOver2 - len2)
            {
              case 0: _score = _mm_extract_epi16( _m2, 0 ); break;
              case 1: _score = _mm_extract_epi16( _m2, 1 ); break;
              case 2: _score = _mm_extract_epi16( _m2, 2 ); break;
              case 3: _score = _mm_extract_epi16( _m2, 3 ); break;
              case 4: _score = _mm_extract_epi16( _m2, 4 ); break;
              case 5: _score = _mm_extract_epi16( _m2, 5 ); break;
              case 6: _score = _mm_extract_epi16( _m2, 6 ); break;
              case 7: _score = _mm_extract_epi16( _m2, 7 ); break;

              default:
                printf("Something is wrong in alignment switch statement\n");
            }

            if (_score < minscore)
            {
                minscore = _score;
                minscoreidx = s+1;
            }
        }

        _m2 = _mm_add_epi16( _m2, _mm_min_epi16( _mm_andnot_si128( _mm_cmpeq_epi16( _seq2win, _seq1win ), _qual2win ), _seq1nqual ) );

        _d2 = _mm_min_epi16( _mm_add_epi16( _d1, _gap_extend ), _mm_add_epi16( _mm_min_epi16( _m1,_i1 ),  // allow I->D
                    _gap_open ) );

        _i2 = _mm_insert_epi16( _mm_srli_si128( _mm_add_epi16( _mm_min_epi16( _mm_add_epi16( _i1,
                                _gap_extend ),
                            _mm_add_epi16( _m1,
                                _gap_open ) ),
                        _nuc_prior ),
                    2 ),
                pos_inf,
                8-1 );
    }

    return (minscore >> 2);
}

//_________________________________________________________________________________________________
