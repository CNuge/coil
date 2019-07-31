

test_that("A short sequence is framed and translated properly", {

  sequence = 'tatgctagggaccgcagttagtgtgattattcgtgctgagttaggacagccaggatcacttattgggaacgatcaaatttacaatacaattgtgactgctcatgcctttattataattttcttcatggtgatacctatcataatcggaggattcggtaattgactggtaccggtaatactaggagcaccagatatagctttccctcgtatgaacaacataagattttgattactccctccttccttaacccttcttataatcgggatactaacagaaagaggggcaggaacaggatgaactgtataccctcctctctcaagaaatatccctcactcaggagctagagtagacctaacaattttttcactacatttagctggagccaggtcacttcttggggctattaatttcatcacaacaattattaatatacgagcagctagaatatctcttgatcgaattcctttatttg'
  sequence_framed = '--------------------------------tatgctagggaccgcagttagtgtgattattcgtgctgagttaggacagccaggatcacttattgggaacgatcaaatttacaatacaattgtgactgctcatgcctttattataattttcttcatggtgatacctatcataatcggaggattcggtaattgactggtaccggtaatactaggagcaccagatatagctttccctcgtatgaacaacataagattttgattactccctccttccttaacccttcttataatcgggatactaacagaaagaggggcaggaacaggatgaactgtataccctcctctctcaagaaatatccctcactcaggagctagagtagacctaacaattttttcactacatttagctggagccaggtcacttcttggggctattaatttcatcacaacaattattaatatacgagcagctagaatatctcttgatcgaattcctttatttg'

  sequence_AAcensored = "-----------MLGTAVSVIIRAELGQPGSLIGNDQIYNTIVTAHAFI?IFFMV?PI?IGGFGNWLVPV?LGAPD?AFPRMNN??FWLLPPSLTLL?IG?LTE?GAGTGWTVYPPLS?NIPHSGA?VDLTIFSLHLAGA?SLLGAINFITTIIN?RAA??SLDRIPLF-"
  sequence_AA5 = "-----------MLGTAVSVIIRAELGQPGSLIGNDQIYNTIVTAHAFIMIFFMVMPIMIGGFGNWLVPVMLGAPDMAFPRMNNMSFWLLPPSLTLLMIGMLTESGAGTGWTVYPPLSSNIPHSGASVDLTIFSLHLAGASSLLGAINFITTIINMRAASMSLDRIPLF"

  dat = coi5p(sequence )
  expect_equal(dat$raw, sequence)
  expect_identical(dat$name, character(0))

  dat = frame(dat)
  expect_equal(dat$framed, sequence_framed)

  dat = translate(dat)
  expect_equal(dat$aaSeq, sequence_AAcensored)

  dat = indel_check(dat)
  expect_equal(dat$indel_likely, TRUE)
  expect_equal(dat$stop_codons, FALSE)

  dat = translate(dat, trans_table = 5)
  expect_equal(dat$aaSeq, sequence_AA5)

  dat = indel_check(dat)
  expect_equal(dat$indel_likely, FALSE)
  expect_equal(dat$stop_codons, FALSE)

})

