

cdef extern from "string.h":
  ctypedef int size_t
  void *memcpy(void *dst,void *src,size_t len)
  void *memmove(void *dst,void *src,size_t len)
  void *memset(void *b,int c,size_t len)

cdef extern from "stdlib.h":
  void free(void *)
  void *malloc(size_t)
  void *calloc(size_t,size_t)
  void *realloc(void *,size_t)
  int c_abs "abs" (int)
  void qsort(void *base, size_t nmemb, size_t size,
             int (*compar)(void *,void *))

cdef extern from "stdio.h":
  ctypedef struct FILE:
    pass
  FILE *fopen(char *,char *)
  FILE *freopen(char *path, char *mode, FILE *stream)
  int fileno(FILE *stream)
  int dup2(int oldfd, int newfd)
  int fflush(FILE *stream)

  FILE * stderr
  FILE * stdout
  int fclose(FILE *)
  int sscanf(char *str,char *fmt,...)
  int printf(char *str,char *fmt,...)
  int sprintf(char *str,char *fmt,...)
  int fprintf(FILE *ifile,char *fmt,...)
  char *fgets(char *str,int size,FILE *ifile)

cdef extern from "ctype.h":
  int toupper(int c)
  int tolower(int c)
  
cdef extern from "unistd.h":
  char *ttyname(int fd)
  int isatty(int fd)  

cdef extern from "string.h":
  int strcmp(char *s1, char *s2)
  int strncmp(char *s1,char *s2,size_t len)
  char *strcpy(char *dest,char *src)
  char *strncpy(char *dest,char *src, size_t len)
  char *strdup(char *)
  char *strcat(char *,char *)
  size_t strlen(char *s)
  int memcmp( void * s1, void *s2, size_t len )

cdef extern from "razf.h":
  pass

cdef extern from "stdint.h":
  ctypedef int int64_t
  ctypedef int int32_t
  ctypedef int uint32_t
  ctypedef int uint8_t
  ctypedef int uint64_t


cdef extern from "bam.h":

  # IF _IOLIB=2, bamFile = BGZF, see bgzf.h
  # samtools uses KNETFILE, check how this works

  ctypedef struct tamFile:
      pass

  ctypedef struct bamFile:
      pass

  ctypedef struct bam1_core_t:
      int32_t tid 
      int32_t pos
      uint32_t bin
      uint32_t qual
      uint32_t l_qname
      uint32_t flag
      uint32_t n_cigar
      int32_t l_qseq
      int32_t mtid 
      int32_t mpos 
      int32_t isize

  ctypedef struct bam1_t:
    bam1_core_t core
    int l_aux
    int data_len
    int m_data
    uint8_t *data

  ctypedef int (*bam_fetch_f)(bam1_t *b, void *data)

  ctypedef struct bam_header_t:
     int32_t n_targets
     char **target_name
     uint32_t *target_len
     void *hash
     void *rg2lib
     int l_text
     char *text

  ctypedef struct bam_index_t:
      pass

  bamFile razf_dopen(int data_fd, char *mode)

  # removed - macros not found

  # int64_t bam_seek( bamFile fp, uint64_t voffset, int where)
  # int64_t bam_tell( bamFile fp )
  # void bam_destroy1( bam1_t * b) 
  # void bam_init_header_hash(bam_header_t *header)

  bam1_t * bam_dup1( bam1_t *src ) 
  
  bam1_t * bam_copy1(bam1_t *bdst, bam1_t *bsrc)
  bam_index_t *bam_index_load(char *f )

  void bam_index_destroy(bam_index_t *idx)

  int bam_parse_region(bam_header_t *header, char *str, int *ref_id, int *begin, int *end)

  int bam_fetch(bamFile fp, bam_index_t *idx, int tid, int beg, int end, void *data, bam_fetch_f func)

  int bam_read1(bamFile fp, bam1_t *b)

  int bam_write1( bamFile fp, bam1_t *b)

  bam_header_t *bam_header_init()

  int bam_header_write( bamFile fp, bam_header_t *header)

  bam_header_t *bam_header_read( bamFile fp )

  void bam_header_destroy(bam_header_t *header)

  bam1_t * bam_dup1( bam1_t *src ) 
  
  bam1_t * bam_copy1(bam1_t *bdst, bam1_t *bsrc)

  uint8_t *bam_aux_get(bam1_t *b,  char tag[2])

  int bam_aux2i(uint8_t *s)
  float bam_aux2f(uint8_t *s)
  double bam_aux2d(uint8_t *s)
  char bam_aux2A( uint8_t *s)
  char *bam_aux2Z( uint8_t *s)
  
  int bam_reg2bin(uint32_t beg, uint32_t end)

  uint32_t bam_calend(bam1_core_t *c, uint32_t *cigar)

cdef extern from "sam.h":

  ctypedef struct samfile_t_un:
    tamFile tamr
    bamFile bam
    FILE *tamw
    
  ctypedef struct samfile_t:
     int type
     samfile_t_un x
     bam_header_t *header

  samfile_t *samopen( char *fn, char * mode, void *aux)

  void samclose(samfile_t *fp)

  int samread(samfile_t *fp, bam1_t *b)

  int samwrite(samfile_t *fp, bam1_t *b)

cdef extern from "pysam_util.h":

    int pysam_dispatch(int argc, char *argv[] )

    # stand-in functions for samtools macros
    void pysam_bam_destroy1( bam1_t * b) 

    # add *nbytes* into the variable length data of *src* at *pos*
    bam1_t * pysam_bam_update( bam1_t * b, 
                               size_t nbytes_old,
                               size_t nbytes_new,
                               uint8_t * pos )

    # translate char to unsigned char
    unsigned char pysam_translate_sequence( char s )

    # stand-ins for samtools macros
    uint32_t * pysam_bam1_cigar( bam1_t * b)
    char * pysam_bam1_qname( bam1_t * b)
    uint8_t * pysam_bam1_seq( bam1_t * b)
    uint8_t * pysam_bam1_qual( bam1_t * b)
    uint8_t * pysam_bam1_aux( bam1_t * b)

    # iterator implemenation
    ctypedef struct bam_fetch_iterator_t:
        pass
  
    bam_fetch_iterator_t* bam_init_fetch_iterator(bamFile fp, bam_index_t *idx, int tid, int beg, int end)
  
    bam1_t * bam_fetch_iterate(bam_fetch_iterator_t *iter)
  
    void bam_cleanup_fetch_iterator(bam_fetch_iterator_t *iter)

###################################################################################################
cdef class Samfile

cdef class IteratorRow:
    cdef bam_fetch_iterator_t*  bam_iter # iterator state object
    cdef bam1_t* b
    cdef error_msg
    cdef int error_state
    cdef Samfile samfile
    cdef bam1_t* getCurrent( self )
    cdef int cnext(self)

###################################################################################################

cdef class IteratorRowAll:
    cdef bam1_t* b
    cdef samfile_t* fp
    cdef bam1_t* getCurrent(self)
    cdef int cnext(self)

###################################################################################################

cdef class Samfile:
    cdef char* filename
    cdef samfile_t* samfile
    cdef bam_index_t *index
    cdef int isbam
    cdef bam1_t* b
    cdef _isOpen( self )
    cdef _hasIndex( self )
    cpdef _open(self, char* filename, mode=*, Samfile template=*, referencenames=*, referencelengths=*, char* text=*, header=*)
    cdef char* getrname(self, int tid)
    cdef _parseRegion(self, reference=*, start=*, end=*, region=*)
    cpdef IteratorRow fetch(self, char* reference, int start, int end)
    cpdef IteratorRowAll fetchAllReads(self)
    cpdef close(self)
    cdef _buildLine(self, fields, record)
    cdef bam_header_t* _buildHeader(self, new_header)
    cdef bam1_t* getCurrent( self )
    cdef int cnext(self)

###################################################################################################

cdef class AlignedRead:
    cdef bam1_t* _delegate
    cdef int hashValue
    cdef int readEnd
    cdef char* _seq
    cdef char* _qual
    cpdef object qname(AlignedRead self)
    cdef char* fastQName(AlignedRead self)
    cdef char* seq(AlignedRead self)
    cdef char* qual(AlignedRead self)
    cdef dict tags(AlignedRead self)
    cdef int flag(AlignedRead self)
    cdef int rname(AlignedRead self)
    cdef int getCigarLength(AlignedRead self)
    cdef int pos(AlignedRead self)
    cdef int end(AlignedRead self)
    cdef bin(AlignedRead self)
    cdef int getCigarOpCode(AlignedRead self, int index)
    cdef int getCigarOpLength(AlignedRead self, int index)
    cdef int rlen(AlignedRead self)
    cpdef int mapq(AlignedRead self)
    cdef int mrnm(AlignedRead self)
    cdef int mpos(AlignedRead self)
    cdef int isize(AlignedRead self)
    cdef int is_paired(AlignedRead self)
    cdef int is_proper_pair(AlignedRead self)
    cdef int is_unmapped(AlignedRead self)
    cdef int mate_is_unmapped(AlignedRead self)
    cdef int is_reverse(AlignedRead self)
    cdef int mate_is_reverse(AlignedRead self)
    cdef int is_read1(AlignedRead self)
    cdef int is_read2(AlignedRead self)
    cdef int is_secondary(AlignedRead self)
    cdef int is_qcfail(AlignedRead self)
    cdef int is_duplicate(AlignedRead self)
    cdef opt(AlignedRead self, tag)

###################################################################################################

cdef AlignedRead makeAlignedRead(bam1_t* src)

###################################################################################################

ctypedef struct cAlignedRead:
    char* seq
    char* qual
    int* cigarOps
    int* cigarLens
    int cigarLen
    int pos
    int rlen
    int end
    int mapq
    int isReverse
    int isPaired
    int isDuplicate
    int isUnMapped
    int insertSize
    int mateIsUnmapped
    int mateIsReverse
    int matePos
    int isReadOne
    int nQualsBelow5
    int nQualsBelow10
    int nQualsBelow15
    int nQualsBelow20

###################################################################################################

cdef cAlignedRead* createRead(bam1_t * src)
cdef destroyRead(cAlignedRead* theRead)
