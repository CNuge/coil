


test_that("A long sequence is framed and translated properly", {

  folmer_seq = ''
  framed_seq = ''
  dat = coi5p_pipe()

  sequence = 'aaccgctgattattttcaaccaaccacaaagatatcggcaaactttatattttatttttggagcttgagctggaatagttggaacatctttaagaattttaattcgagctgaattaggacatcctggagcattaattggagatgatcaaatttataatgtaattgtaactgcacatgcttttattataattttttttatggttatacctattataattggtggatttggaaattgattagtgcctttaatattaggtgctcctgatatagcattcccacgaataaataatataagattttgactactacctcctgctctttctttactattagtaagtagaatagttgaaaatggagctggaacaggatgaactgtttatccacctttatccgctggaattgctcatggtggagcttcagttgatttagctattttttctctacatttagcagggatttcttcaattttaggagctctaaattttattacaactgtaattaatatacgatcaacaggaatttcattagatcgtatacctttatttgtttgatcagtagttattactgctttattattgttattatcacttccagtactagcaggagctattactatattattaacagatcgaaatttaaatacatcattttttgacccagcgggaggaggagatcctattttatatcaacatttattatttattta'
  sequence_framed = '-ctttatattttatttttggagcttgagctggaatagttggaacatctttaagaattttaattcgagctgaattaggacatcctggagcattaattggagatgatcaaatttataatgtaattgtaactgcacatgcttttattataattttttttatggttatacctattataattggtggatttggaaattgattagtgcctttaatattaggtgctcctgatatagcattcccacgaataaataatataagattttgactactacctcctgctctttctttactattagtaagtagaatagttgaaaatggagctggaacaggatgaactgtttatccacctttatccgctggaattgctcatggtggagcttcagttgatttagctattttttctctacatttagcagggatttcttcaattttaggagctctaaattttattacaactgtaattaatatacgatcaacaggaatttcattagatcgtatacctttatttgtttgatcagtagttattactgctttattattgttattatcacttccagtactagcaggagctattactatattattaacagatcgaaatttaaatacatcattttttgacccagcgggaggaggagatcctattttatatcaacattta'

  sequence_AAcensored = "-LYFIFGAWAG?VGTSL?ILIRAELGHPGALIGDDQIYNVIVTAHAFI?IFFMV?PI?IGGFGNWLVPL?LGAPD?AFPR?NN??FWLLPPALSLLLVS??VENGAGTGWTVYPPLSAGIAHGGASVDLAIFSLHLAGISSILGALNFITTVIN?RSTGISLDR?PLFVWSVVITALLLLLSLPVLAGAIT?LLTDRNLNTSFFDPAGGGDPILYQHL"
  sequence_AA5 = "-LYFIFGAWAGMVGTSLSILIRAELGHPGALIGDDQIYNVIVTAHAFIMIFFMVMPIMIGGFGNWLVPLMLGAPDMAFPRMNNMSFWLLPPALSLLLVSSMVENGAGTGWTVYPPLSAGIAHGGASVDLAIFSLHLAGISSILGALNFITTVINMRSTGISLDRMPLFVWSVVITALLLLLSLPVLAGAITMLLTDRNLNTSFFDPAGGGDPILYQHL"

  dat = coi5p(sequence)
  expect_equal(dat$raw, sequence)
  expect_identical(dat$name, character(0))

  dat = frame(dat)
  expect_equal(dat$framed, sequence_framed)

  dat = translate(dat)
  expect_equal(dat$aaSeq, sequence_AAcensored)

  dat = indel_check(dat)
  expect_equal(dat$indel_likely, FALSE)
  expect_equal(dat$stop_codons, FALSE)

  dat = translate(dat, trans_table = 5)
  expect_equal(dat$aaSeq, sequence_AA5)

  dat = indel_check(dat)
  expect_equal(dat$indel_likely, FALSE)
  expect_equal(dat$stop_codons, FALSE)

})
