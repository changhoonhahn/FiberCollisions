import numpy as np
import sys
''' 
Simple python script that combines 0 to n (input) pthalo dr11 v11p0 random files
'''
if __name__=="__main__":
    pthalo_dir="/mount/riachuelo1/hahn/data/PTHalo/v11p0/"
    pthalo_file_prefix="cmass_dr11_north_randoms_ir4"
    pthalo_file_suffix=".v11.0.wghtv.txt"
    n = sys.argv[1]
    print ''.join(["COMBINING ", str(n), " RANDOM"])

    pthalo_randcombine_filename = ''.join([pthalo_dir, 'cmass_dr11_north_', str(n), '_randoms_ir4_combined.v11.0.txt']) 
    pthalo_randcombine_file = open(pthalo_randcombine_filename,'w')

    for i in range(1,int(n)+1):
        pthalo_rand_filename = ''.join([pthalo_dir, pthalo_file_prefix, str(i+1000)[1:4], pthalo_file_suffix])
        data = np.loadtxt(pthalo_rand_filename)

        for j in range(len(data)):
            pthalo_randcombine_file.write(
                    str(data[j,0])+'\t'+str(data[j,1])+'\t'+str(data[j,2])+'\t'+str(data[j,3])+'\t'+
                    str(data[j,4])+'\t'+str(data[j,5])+'\t'+str(data[j,6])+'\t'+str(data[j,7])+'\n')
