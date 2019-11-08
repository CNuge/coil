test_that("The reading frame of sequences is reported properly", {

  #TODO - add $align_report description to the readme and the vignette

  long_sequence = 'aaccgctgattattttcaaccaaccacaaagatatcggcaaactttatattttatttttggagcttgagctggaatagttggaacatctttaagaattttaattcgagctgaattaggacatcctggagcattaattggagatgatcaaatttataatgtaattgtaactgcacatgcttttattataattttttttatggttatacctattataattggtggatttggaaattgattagtgcctttaatattaggtgctcctgatatagcattcccacgaataaataatataagattttgactactacctcctgctctttctttactattagtaagtagaatagttgaaaatggagctggaacaggatgaactgtttatccacctttatccgctggaattgctcatggtggagcttcagttgatttagctattttttctctacatttagcagggatttcttcaattttaggagctctaaattttattacaactgtaattaatatacgatcaacaggaatttcattagatcgtatacctttatttgtttgatcagtagttattactgctttattattgttattatcacttccagtactagcaggagctattactatattattaacagatcgaaatttaaatacatcattttttgacccagcgggaggaggagatcctattttatatcaacatttattatttattta'
  long_sequence_framed = '-ctttatattttatttttggagcttgagctggaatagttggaacatctttaagaattttaattcgagctgaattaggacatcctggagcattaattggagatgatcaaatttataatgtaattgtaactgcacatgcttttattataattttttttatggttatacctattataattggtggatttggaaattgattagtgcctttaatattaggtgctcctgatatagcattcccacgaataaataatataagattttgactactacctcctgctctttctttactattagtaagtagaatagttgaaaatggagctggaacaggatgaactgtttatccacctttatccgctggaattgctcatggtggagcttcagttgatttagctattttttctctacatttagcagggatttcttcaattttaggagctctaaattttattacaactgtaattaatatacgatcaacaggaatttcattagatcgtatacctttatttgtttgatcagtagttattactgctttattattgttattatcacttccagtactagcaggagctattactatattattaacagatcgaaatttaaatacatcattttttgacccagcgggaggaggagatcctattttatatcaacattta'

  #expected:
  #42 leading bp trimmed and there is a leading dash
  x_long = 43
  y_long = 2

  long_desc = paste0("bp ", x_long, " of the raw sequence is base pair ", y_long, " of the Folmer region.")

  long_dat = coi5p(long_sequence, name = 'long_frame_report')
  long_dat = frame(long_dat)

  long_dat$data$raw_int_trim
  long_dat$data$raw_start
  long_dat$data$folmer_start
  long_dat$align_report
  long_dat$was_trimmed

  # long_desc == long_dat$align_report
  expect_equal(long_desc, long_dat$align_report)

  short_sequence = 'tatgctagggaccgcagttagtgtgattattcgtgctgagttaggacagccaggatcacttattgggaacgatcaaatttacaatacaattgtgactgctcatgcctttattataattttcttcatggtgatacctatcataatcggaggattcggtaattgactggtaccggtaatactaggagcaccagatatagctttccctcgtatgaacaacataagattttgattactccctccttccttaacccttcttataatcgggatactaacagaaagaggggcaggaacaggatgaactgtataccctcctctctcaagaaatatccctcactcaggagctagagtagacctaacaattttttcactacatttagctggagccaggtcacttcttggggctattaatttcatcacaacaattattaatatacgagcagctagaatatctcttgatcgaattcctttatttg'
  short_sequence_framed = '--------------------------------tatgctagggaccgcagttagtgtgattattcgtgctgagttaggacagccaggatcacttattgggaacgatcaaatttacaatacaattgtgactgctcatgcctttattataattttcttcatggtgatacctatcataatcggaggattcggtaattgactggtaccggtaatactaggagcaccagatatagctttccctcgtatgaacaacataagattttgattactccctccttccttaacccttcttataatcgggatactaacagaaagaggggcaggaacaggatgaactgtataccctcctctctcaagaaatatccctcactcaggagctagagtagacctaacaattttttcactacatttagctggagccaggtcacttcttggggctattaatttcatcacaacaattattaatatacgagcagctagaatatctcttgatcgaattcctttattt'

  x_short = 1
  y_short = 33
  short_desc = paste0("bp ", x_short, " of the raw sequence is base pair ", y_short, " of the Folmer region.")

  short_dat = coi5p(short_sequence, name = 'short_frame_report')
  short_dat = frame(short_dat)

  short_dat$data$raw_start
  short_dat$data$folmer_start

  # short_desc == short_dat$align_report
  expect_equal(short_desc, short_dat$align_report)


  normal_sequence = 'ctttacctgatttttggtgcatgagcaggtatagttggaacagccctaagtctcctaattcgagctgaacttgggcaacctggatcacttttaggagatgatcagatttataatgtaatcgtaaccgcccacgcttttgtaataatctttttcatggttataccaattataattggtggtttcggaaattgattagttcctttaataattggagcgccagatatagccttcccacgaataaataacataagtttctgacttcttccaccatcatttcttcttctcctcgcctctgctggagtagaagctggagcaggtactggttgaacagtttatcctccattagctagcaatctagcacatgctggaccatctgttgatttagctattttttctcttcacttagccggtgtttcatcaattttagcttcaattaattttatcacaaccattattaatataaaaccaccagctatttcccaatatcaaacaccattatttgtttgatctattcttgtaaccactattcttcttctcctctcacttccagttcttgcagcaggaattacaatattacttacagatcgtaaccttaatactacattctttgaccctgcaggtggaggagacccaatcctttatcaacattta'
  normal_sequence_framed = '---ctttacctgatttttggtgcatgagcaggtatagttggaacagccctaagtctcctaattcgagctgaacttgggcaacctggatcacttttaggagatgatcagatttataatgtaatcgtaaccgcccacgcttttgtaataatctttttcatggttataccaattataattggtggtttcggaaattgattagttcctttaataattggagcgccagatatagccttcccacgaataaataacataagtttctgacttcttccaccatcatttcttcttctcctcgcctctgctggagtagaagctggagcaggtactggttgaacagtttatcctccattagctagcaatctagcacatgctggaccatctgttgatttagctattttttctcttcacttagccggtgtttcatcaattttagcttcaattaattttatcacaaccattattaatataaaaccaccagctatttcccaatatcaaacaccattatttgtttgatctattcttgtaaccactattcttcttctcctctcacttccagttcttgcagcaggaattacaatattacttacagatcgtaaccttaatactacattctttgaccctgcaggtggaggagacccaatcctttatcaacattta'

  x_norm = 1
  y_norm = 4
  normal_desc = paste0("bp ", x_norm, " of the raw sequence is base pair ", y_norm, " of the Folmer region.")

  normal_dat = coi5p(normal_sequence, name = 'normal_frame_report')
  normal_dat = frame(normal_dat)

  # normal_desc == normal_dat$align_report
  expect_equal(normal_desc, normal_dat$align_report)


  #  COI is defined relative to the mouse mitochondrial genome as the
  #  648 bp region that starts at position 58 and stops at position 705
  #  The PHMM only spans the first 657bp of the sequence

})
