import ftputil
import os

filenames = [filename for filename in os.listdir('.') if filename.endswith('fastq.gz')]

with ftputil.FTPHost('ftp-private.ncbi.nlm.nih.gov', 'subftp', 'SniappegEtnurak3') as host:
    print(host.getcwd())
    host.chdir('uploads/dwc001_ucsd.edu_bdLPdXnY/Methyldeoxaphomin_NPDG_F-Dd2pol-resistant_clones_2')
    print(host.getcwd())
    for filename in filenames:
        print("Uploading %s" % filename)
        host.upload(filename, filename)
    print("----")
    print("Done")
