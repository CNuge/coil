

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

  #####
  #the following are additional tests that must be passed to ensure coi5p is working well
  #coi5p barcode sequence: CFF077-16
  sequence_2 = 'attatacctcctatttggggcctgagcaggaatagtagggacagccctaagcatcctaattcgagcagaacttggacagccaggcgctctattaggtgacgaccaaatttataacgtagtcgttacagcccacgcattcgtcataattttcttcatagtaataccaataataatcggcgggttcggtaactgattagtcccactaataattggagcccccgacatagcattcccacgaataaacaatataagcttctgactcctgcccccatcattccttctacttttagcctcatctatagtagaagcaggggccggaacaggatgaactgtttacccaccactagccggcaatttagcacacgcaggagcatccgtagacctaaccattttttccctacacctagcaggagtttcatcaattttaggcgctatcaatttcattactaccattatcaatataaaaccgccagccatatcacaatatcaaacccccctatttgtgtgatcagtccttattaccgccgtactacttctcctctccctcccagtcttagccgcaggcattacaatactcttaaccgaccgaaatctaaatactacattcttcgacccagccggaggcggtgaccccattttataccaacatcta'
  #trimmed version of this, two missing from front
           # "acattatacctcctatttggggcctgagcaggaatagtagggacagccctaagcatcctaattcgagcagaacttggacagccaggcgctctattaggtgacgaccaaatttataacgtagtcgttacagcccacgcattcgtcataattttcttcatagtaataccaataataatcggcgggttcggtaactgattagtcccactaataattggagcccccgacatagcattcccacgaataaacaatataagcttctgactcctgcccccatcattccttctacttttagcctcatctatagtagaagcaggggccggaacaggatgaactgtttacccaccactagccggcaatttagcacacgcaggagcatccgtagacctaaccattttttccctacacctagcaggagtttcatcaattttaggcgctatcaatttcattactaccattatcaatataaaaccgccagccatatcacaatatcaaacccccctatttgtgtgatcagtccttattaccgccgtactacttctcctctccctcccagtcttagccgcaggcattacaatactcttaaccgaccgaaatctaaatactacattcttcgacccagccggaggcggtgaccccattttataccaacatcta"
  sequence2_framed = ''

  sequence2_AAcensored = ''
  sequence2_AA5 = ''


  dat2 = coi5p(sequence )
  expect_equal(dat2$raw, sequence)
  expect_identical(dat2$name, character(0))

  dat2 = frame(dat2)
  expect_equal(dat2$framed, sequence_framed)

  dat2 = translate(dat2)
  expect_equal(dat2$aaSeq, sequence_AAcensored)

  dat2 = indel_check(dat2)
  expect_equal(dat2$indel_likely, TRUE)
  expect_equal(dat2$stop_codons, FALSE)

  dat2 = translate(dat2, trans_table = 2)
  expect_equal(dat2$aaSeq, sequence_AA5)

  dat2 = indel_check(dat2)
  expect_equal(dat2$indel_likely, FALSE)
  expect_equal(dat2$stop_codons, FALSE)






  sequence_3 = 'ggcgctcttctgggggatgaccaaatctataacgtgatcgtcacagcccatgccttcgttatgattttctttatagtcatgccaattataatcgggggctttggaaactgattaattcccctaataatcggagcccctgatatggcattccctcgaataaataacataagcttctgactccttcctccatcctttctcctcctcctgtcttcatcaggagttgaagccggcgcgggtactggatgaacagtatacccccctctagccggcaacctcgcccacgcaggagcctctgttgatttaactatcttctcccttcatttagctggaatctcctcaattttaggagccattaattttattacgaccattattaacataaaacctccagccatctctcagtaccaaaccccccttttcgtttgagccgtgctagttactgctgtccttctattactttccctccccgtcctggca'
  #trimmed version of this, 81 missing from front, 93 missing from back.
  #"ctctatttagtatttggtgcctgagccgggatagtaggcaccgccctgagtctactgattcgggcggaactaagccagccgggcgctcttctgggggatgaccaaatctataacgtgatcgtcacagcccatgccttcgttatgattttctttatagtcatgccaattataatcgggggctttggaaactgattaattcccctaataatcggagcccctgatatggcattccctcgaataaataacataagcttctgactccttcctccatcctttctcctcctcctgtcttcatcaggagttgaagccggcgcgggtactggatgaacagtatacccccctctagccggcaacctcgcccacgcaggagcctctgttgatttaactatcttctcccttcatttagctggaatctcctcaattttaggagccattaattttattacgaccattattaacataaaacctccagccatctctcagtaccaaaccccccttttcgtttgagccgtgctagttactgctgtccttctattactttccctccccgtcctggcagcaggcattactatgttacttacagaccgaaatctaaacaccactttctttgacccggcaggcgggggagatccaattttataccaacacctc"
  sequence3_framed = ''

  sequence3_AAcensored = ''
  sequence3_AA5 = ''


  dat3 = coi5p(sequence )
  expect_equal(dat3$raw, sequence)
  expect_identical(dat3$name, character(0))

  dat3 = frame(dat3)
  expect_equal(dat3$framed, sequence_framed)

  dat3 = translate(dat3)
  expect_equal(dat3$aaSeq, sequence_AAcensored)

  dat3 = indel_check(dat3)
  expect_equal(dat3$indel_likely, TRUE)
  expect_equal(dat3$stop_codons, FALSE)

  dat3 = translate(dat3, trans_table = 2)
  expect_equal(dat3$aaSeq, sequence_AA5)

  dat3 = indel_check(dat3)
  expect_equal(dat3$indel_likely, FALSE)
  expect_equal(dat3$stop_codons, FALSE)


})

