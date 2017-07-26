from libc.stdlib cimport malloc, free


cdef int MINIMUM_TANDEM_LENGTH = 4
cdef dict repunit_dict = None


cdef extern from "tandem.h":
    void annotate( char* sequence, char* sizes, char* displacements, int length)


cdef tuple calculate_size_and_displacement( bytes sequence ):
    cdef char* csequence = sequence
    cdef int length = len(sequence)
    cdef char* csizes = <char*>malloc(length+1)
    cdef char* cdisplacements = <char*>malloc(length+1)
    annotate( csequence, csizes, cdisplacements, length )
    cdef bytes sizes = csizes
    cdef bytes displacements = cdisplacements
    return (sizes, displacements)


cdef bytes normalize_repunit_fast( bytes unit ):
    global repunit_dict
    cdef int length = len(unit)
    cdef char* cunit = unit
    cdef int both
    if length == 1:
        if cunit[0] <= 'C': return unit
        if cunit[0] == 'G': return bytes('c')
        return bytes('g')
    elif length == 2:
        both = <int>cunit[0]+<int>cunit[1]
        if both == 138: return bytes('CG')   # c+g
        elif both == 149: return bytes('AT') # a+t
        elif both == 132: return bytes('AC') # a+c 
        elif both == 155: return bytes('ac') # g+t
        elif both == 136: return bytes('AG') # a+g
        return bytes('ag') # t+c
    cdef bytes n1
    cdef bytes n2
    cdef bytes n3
    cdef bytes n4
    cdef bytes n5
    cdef bytes n6
    if repunit_dict == None:
        repunit_dict = {}
        for n1 in bytes('ACGT'):
            for n2 in bytes('ACGT'):
                for n3 in bytes('ACGT'):
                    repunit_dict[n1+n2+n3] = normalize_repunit(n1+n2+n3)
                    for n4 in bytes('ACGT'):
                        repunit_dict[n1+n2+n3+n4] = normalize_repunit(n1+n2+n3+n4)
                        for n5 in bytes('ACGT'):
                            repunit_dict[n1+n2+n3+n4+n5] = normalize_repunit(n1+n2+n3+n4+n5)
                            for n6 in bytes('ACGT'):
                                repunit_dict[n1+n2+n3+n4+n5+n6] = normalize_repunit(n1+n2+n3+n4+n5+n6)
    if length <= 6:
        return repunit_dict[unit]
    return normalize_repunit(unit)


cdef bytes normalize_repunit( bytes unit ):
    cdef list unit2list = [{'A':'T','T':'A','C':'G','G':'C'}.get(c,'N') for c in unit]
    unit2list.reverse()
    cdef bytes unit2 = bytes(''.join(bytes(x) for x in unit2list))
    unit2 += unit2
    unit += unit
    cdef int length = len(unit)
    cdef bytes normunit = sorted([ unit[i:i+length] for i in range(length) ] + [ unit2[i:i+length]+'-' for i in range(length) ])[0]
    if normunit[-1] == '-':
        normunit = normunit[:-1].lower()
    return normunit


cdef tuple try_microsat( char* seq, int seqlen, int start, int replen, int direction ):
    cdef int minpos = start
    cdef int maxpos = start
    cdef int length
    cdef bytes character
    while minpos >= 0 and 0 <= minpos+replen*direction < seqlen and seq[minpos] == seq[minpos + replen*direction] and seq[minpos] != 'N': 
        minpos -= 1
    while maxpos < seqlen and maxpos+replen*direction < seqlen and seq[maxpos] == seq[maxpos + replen*direction] and seq[maxpos] != 'N': 
        maxpos += 1
    length = maxpos - minpos + replen - 1
    if length < 2*replen:
        character = seq[start]
        return 1, character, start
    return length, seq[minpos+1:minpos+replen+1], min(minpos+1,minpos+1+replen*direction)


cdef microsattelite( char* seq, int seqlen, int start, int tandem ):
    cdef tuple maxlen = (-1,"",-1)
    cdef tuple result
    for m in range(1,tandem+1):
        result = try_microsat(seq, seqlen, start, m, -1)

        if result[0] > maxlen[0]: 
            maxlen = result
        result = try_microsat(seq, seqlen, start, m, 1)

        if result[0] > maxlen[0]: 
            maxlen = result
    return maxlen


def get_annotation(seq, pos, tandem):
    t, repunitT, repstart = microsattelite( seq, len(seq), pos, tandem )
    return t,normalize_repunit(repunitT),repstart


def add_tandem(pos, tandemlen, tandemunit, indelq, indel_q_data, output_base=0):
    """ Adds gap opening penalties in string indelq, for the tandem repeat described by output.
        indel_q_data contains a dictionary of error models.  These have the following three forms:
        {'nonrepetitive' : '#',   # must be a single character
         'AAC' : '#######',       # gap opening penalties, phred 33-based, for AAC repeats of length 1,2,..,7
         3 : '######'}            # gap opening penalties for any other triplet repeats of length 1,2,...,6 """
    tandemunit = tandemunit.upper()
    if pos == -1: 
        return
    if tandemunit in indel_q_data:
        qdata = indel_q_data[tandemunit]
    elif len(tandemunit) in indel_q_data:
        qdata = indel_q_data[len(tandemunit)]
    else:
        return
    if tandemlen-1 >= len(qdata): 
        q = chr(ord(qdata[-1])-ord('!')+output_base)
    else:
        q = chr(ord(qdata[tandemlen-1])-ord('!')+output_base)
    for idx in range(pos, pos+tandemlen):
        indelq[idx] = min(q, indelq[idx])


cdef bytes annotate_sequence_slow(sequence, dict indel_q_data, int output_base, int tandem):
    """ Annotates a sequence with the local indel probability (gap opening penalty),
        using an error model described in indel_q_data """
    pos = 0
    indelq = [chr(ord(indel_q_data['nonrepetitive'])-ord('!')+output_base)] * len(sequence)
    oldoutput = [-1,-1,-1]
    
    while pos < len(sequence):
        tandemlen, tandemunit, repstart = get_annotation(sequence, pos, tandem)

        output = [repstart, tandemlen, tandemunit]

        # concatenate or output
        if tandemlen >= 2:
            if oldoutput[0] + oldoutput[1] >= output[0] and oldoutput[2] == output[2]:
                oldoutput[1] = output[0] + output[1] - oldoutput[0]
            else:
                if oldoutput[0] != -1:
                    add_tandem(oldoutput[0], oldoutput[1], oldoutput[2], indelq, indel_q_data, output_base = output_base)
                oldoutput = output
            pos = repstart + tandemlen
        else:
            pos += 1
    if oldoutput[0] != -1:
        add_tandem(oldoutput[0], oldoutput[1], oldoutput[2], indelq, indel_q_data, output_base = output_base)
    return bytes(''.join(bytes(x) for x in indelq))


cdef bytes annotate_sequence(sequence, dict indel_q_data, int output_base, int tandem):
    """ Annotates a sequence with the local indel probability (gap opening penalty),
        using an error model described in indel_q_data """

    cdef list indelq = [chr(ord(indel_q_data['nonrepetitive'])-ord('!')+output_base)] * len(sequence)

    cdef bytes sizes
    cdef bytes displacements
    (sizes, displacements) = calculate_size_and_displacement( sequence )

    cdef char* csizes = sizes
    cdef char* cdisplacements = displacements
    cdef int length = len(sequence)
    cdef int pos = 0
    cdef int tandemunitlen
    cdef int tandemlen
    cdef bytes tandemunit

    cdef int oldoutputpos = -1
    cdef int oldoutputlen = -1
    cdef bytes oldoutputunit

    while pos < length:

        tandemunitlen = cdisplacements[pos]
        tandemlen = csizes[pos]
        tandemunit = normalize_repunit_fast( sequence[pos:pos+tandemunitlen] )

        # concatenate or output
        if tandemlen >= 2:
            if oldoutputpos + oldoutputlen >= pos and oldoutputunit == tandemunit:
                oldoutputlen = pos + tandemlen - oldoutputpos
            else:
                if oldoutputpos != -1 and oldoutputlen >= MINIMUM_TANDEM_LENGTH:
                    add_tandem( oldoutputpos, oldoutputlen, oldoutputunit, indelq, indel_q_data, output_base = output_base)
                oldoutputpos = pos
                oldoutputlen = tandemlen
                oldoutputunit = tandemunit
        
        # next position
        pos += 1

    # add last tandem
    if oldoutputpos != -1:
        add_tandem( oldoutputpos, oldoutputlen, oldoutputunit, indelq, indel_q_data, output_base = output_base)

    # return result
    return bytes(''.join(bytes(x) for x in indelq))


cpdef testAnnotate():
    seq1 = "TATTTGCATGCGCTTTCGAGCTGTTGAAGAGACGTGTATTGGAATAAGTAATCACATAAGTGTTAGTAACTTATTTAAATACGTATAGAGTCGCCTATTTGCCTAGCCTTTTGGTTCTCAGATTTTTTAATTATTACATTGCTATAAGGGTGTAACTGTGTGATAGCCAAAATTTTAAGCTGCAAATGGTTTGTAAATATGATATATTACAAGCTTCATGAAAATCGGTTTATGACTGATCCGCGATTACGTTGAAAGGCGACTGGCAGAGATACTTTTGTTCAGATGTTTTTTCAGGTAGCGATTCCAATGAATAGGTAAAATACCTTGCAAGTTTTGTTGTTGTCGTTGGAGGAAATGTGGATGTGGTTGTTATTGTTGA"                      
    
    import time

    t = time.clock()
    model = {'nonrepetitive':'S',
             1:'SSI?5+#',
             'AG':'SS#'}
    indelq = annotate_sequence_slow(seq1, model, ord('!'), 24)
    indelq_fast = annotate_sequence(seq1, model, ord('!'), 24)
    print seq1
    print indelq
    print indelq_fast

    # test speed
    seq2 = seq1*10
    for i in range(10000):
        indelq = annotate_sequence(seq1, model, ord('!'), 24)

    print time.clock()-t


