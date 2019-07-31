
test_that("A normal sequence in framed and translated properly", {

  sequence = 'ctttacctgatttttggtgcatgagcaggtatagttggaacagccctaagtctcctaattcgagctgaacttgggcaacctggatcacttttaggagatgatcagatttataatgtaatcgtaaccgcccacgcttttgtaataatctttttcatggttataccaattataattggtggtttcggaaattgattagttcctttaataattggagcgccagatatagccttcccacgaataaataacataagtttctgacttcttccaccatcatttcttcttctcctcgcctctgctggagtagaagctggagcaggtactggttgaacagtttatcctccattagctagcaatctagcacatgctggaccatctgttgatttagctattttttctcttcacttagccggtgtttcatcaattttagcttcaattaattttatcacaaccattattaatataaaaccaccagctatttcccaatatcaaacaccattatttgtttgatctattcttgtaaccactattcttcttctcctctcacttccagttcttgcagcaggaattacaatattacttacagatcgtaaccttaatactacattctttgaccctgcaggtggaggagacccaatcctttatcaacattta'
  sequence_framed = '---ctttacctgatttttggtgcatgagcaggtatagttggaacagccctaagtctcctaattcgagctgaacttgggcaacctggatcacttttaggagatgatcagatttataatgtaatcgtaaccgcccacgcttttgtaataatctttttcatggttataccaattataattggtggtttcggaaattgattagttcctttaataattggagcgccagatatagccttcccacgaataaataacataagtttctgacttcttccaccatcatttcttcttctcctcgcctctgctggagtagaagctggagcaggtactggttgaacagtttatcctccattagctagcaatctagcacatgctggaccatctgttgatttagctattttttctcttcacttagccggtgtttcatcaattttagcttcaattaattttatcacaaccattattaatataaaaccaccagctatttcccaatatcaaacaccattatttgtttgatctattcttgtaaccactattcttcttctcctctcacttccagttcttgcagcaggaattacaatattacttacagatcgtaaccttaatactacattctttgaccctgcaggtggaggagacccaatcctttatcaacattta'
  seqname = 'test_seq1'
  sequence_AAcensored = "-LYLIFGAWAG?VGTALSLLIRAELGQPGSLLGDDQIYNVIVTAHAFV?IFFMV?PI?IGGFGNWLVPL?IGAPD?AFPR?NN?SFWLLPPSFLLLLASAGVEAGAGTGWTVYPPLASNLAHAGPSVDLAIFSLHLAGVSSILASINFITTIIN??PPAISQYQTPLFVWSILVTTILLLLSLPVLAAGIT?LLTDRNLNTTFFDPAGGGDPILYQHL"
  sequence_AA5 ="-LYLIFGAWAGMVGTALSLLIRAELGQPGSLLGDDQIYNVIVTAHAFVMIFFMVMPIMIGGFGNWLVPLMIGAPDMAFPRMNNMSFWLLPPSFLLLLASAGVEAGAGTGWTVYPPLASNLAHAGPSVDLAIFSLHLAGVSSILASINFITTIINMKPPAISQYQTPLFVWSILVTTILLLLSLPVLAAGITMLLTDRNLNTTFFDPAGGGDPILYQHL"
  dat = coi5p(sequence, name = seqname)
  expect_equal(dat$raw, sequence)
  expect_equal(dat$name, seqname)


  dat = frame(dat)
  expect_equal(dat$framed, sequence_framed)

  dat = translate(dat)
  expect_equal(dat$aaSeq, sequence_AAcensored)

  dat = translate(dat, trans_table = 5)
  expect_equal(dat$aaSeq, sequence_AA5)

  dat = indel_check(dat)
  expect_equal(dat$indel_likely, FALSE)
  expect_equal(dat$stop_codons, FALSE)

})



context("A normal coi5p sequence")

seq_normal =


test_normal = coi5p(seq_normal)
test_normal
test_normal = frame(test_normal)
test_normal
test_normal = translate(test_normal)
test_normal
test_normal = translate(test_normal, trans_table = 5)
test_normal
test_normal = indel_check(test_normal)
test_normal

#once the tests are written above, copy them into the other files and interface with
#sequence data that matches the different criteria given. don't generate the
#expected using coi5p library (back to scripts/functions/raw data) so that the tests aren't tautological.
