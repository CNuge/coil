test_that("coi5p objects are printed properly", {

  #baseline obj
  seqname = "test_seq1"
  sequence = 'ctttacctgatttttggtgcatgagcaggtatagttggaacagccctaagtctcctaattcgagctgaacttgggcaacctggatcacttttaggagatgatcagatttataatgtaatcgtaaccgcccacgcttttgtaataatctttttcatggttataccaattataattggtggtttcggaaattgattagttcctttaataattggagcgccagatatagccttcccacgaataaataacataagtttctgacttcttccaccatcatttcttcttctcctcgcctctgctggagtagaagctggagcaggtactggttgaacagtttatcctccattagctagcaatctagcacatgctggaccatctgttgatttagctattttttctcttcacttagccggtgtttcatcaattttagcttcaattaattttatcacaaccattattaatataaaaccaccagctatttcccaatatcaaacaccattatttgtttgatctattcttgtaaccactattcttcttctcctctcacttccagttcttgcagcaggaattacaatattacttacagatcgtaaccttaatactacattctttgaccctgcaggtggaggagacccaatcctttatcaacattta'

  dat = coi5p(sequence)

  expected1 = "coi5p barcode sequence
raw sequence:
ctttacctgatttttggtgcatgag...agacccaatcctttatcaacattta"

  dstr1 = capture_output(dat, print=TRUE)

  expect_equal(dstr1, expected1)


  expected2 = "coi5p barcode sequence: test_seq1
raw sequence:
ctttacctgatttttggtgcatgag...agacccaatcctttatcaacattta
framed sequence:
---ctttacctgatttttggtgcat...agacccaatcctttatcaacattta
Amino acid sequence:
-LYLIFGAWAG?VGTALSLLIRAEL...LTDRNLNTTFFDPAGGGDPILYQHL
The sequence likely does not contain an insertion or deletion.
Stop codon present: FALSE, Amino acid PHMM score:-198.41257"

#obj with name

dat = coi5p_pipe(sequence, name = seqname)

expected2 = "coi5p barcode sequence: test_seq1
raw sequence:
ctttacctgatttttggtgcatgag...agacccaatcctttatcaacattta
framed sequence:
---ctttacctgatttttggtgcat...agacccaatcctttatcaacattta
Amino acid sequence:
-LYLIFGAWAG?VGTALSLLIRAEL...LTDRNLNTTFFDPAGGGDPILYQHL
The sequence likely does not contain an insertion or deletion.
Stop codon present: FALSE, Amino acid PHMM score:-198.41257"

dstr2 = capture_output(dat, print=TRUE)

expect_equal(dstr2, expected2)


  #object after each stage
  dat = coi5p(sequence, name = seqname)

  expected3 = "coi5p barcode sequence: test_seq1
raw sequence:
ctttacctgatttttggtgcatgag...agacccaatcctttatcaacattta"

  dstr3 = capture_output(dat, print=TRUE)

  expect_equal(dstr3, expected3)



  dat = frame(dat)
  expected4 = "coi5p barcode sequence: test_seq1
raw sequence:
ctttacctgatttttggtgcatgag...agacccaatcctttatcaacattta
framed sequence:
---ctttacctgatttttggtgcat...agacccaatcctttatcaacattta"

  dstr4 = capture_output(dat, print=TRUE)

  expect_equal(dstr4, expected4)




  dat = translate(dat)

  expected5 = "coi5p barcode sequence: test_seq1
raw sequence:
ctttacctgatttttggtgcatgag...agacccaatcctttatcaacattta
framed sequence:
---ctttacctgatttttggtgcat...agacccaatcctttatcaacattta
Amino acid sequence:
-LYLIFGAWAG?VGTALSLLIRAEL...LTDRNLNTTFFDPAGGGDPILYQHL"

  dstr5 = capture_output(dat, print=TRUE)

  expect_equal(dstr5, expected5)


  dat = indel_check(dat)

  expected6 = "coi5p barcode sequence: test_seq1
raw sequence:
ctttacctgatttttggtgcatgag...agacccaatcctttatcaacattta
framed sequence:
---ctttacctgatttttggtgcat...agacccaatcctttatcaacattta
Amino acid sequence:
-LYLIFGAWAG?VGTALSLLIRAEL...LTDRNLNTTFFDPAGGGDPILYQHL
The sequence likely does not contain an insertion or deletion.
Stop codon present: FALSE, Amino acid PHMM score:-198.41257"

  dstr6 = capture_output(dat, print=TRUE)

  expect_equal(dstr6, expected6)

})
