import lbcs_stats

#picdir = '/data/data-disk1/data/LBCS/lbcs_data/picfiles/' ### pic dir
picdir = '/data/data-disk1/data/LBCS/lbcs_testdata/' ### pic dir

#lbcs_stats.getstats(picdir) ### generate stats
lbcs_stats.plotstats() ### generates plot
lbcs_stats.plotit()
lbcs_stats.plotdetect()

print '\n\n*** Stats Ready! ***\n'




