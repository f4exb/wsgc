echo "*** Message correlation OOK" | tee validate.out
/opt/install/wsgc_test/bin/wsgc_test -s 4096 -c 1023 -C 10 -t 0.0 -r 0.0 -N 4 -I 0 -p 1,3,5,7,9,11,13,15 -B 2 -z 4 -n -6 -d OOK -H 0.13,1.1 >> validate.out
echo "*** Message correlation DBPSK" | tee -a validate.out
/opt/install/wsgc_test/bin/wsgc_test -s 4096 -c 1023 -C 10 -t 0.0 -r 0.0 -N 4 -I 0 -p 1,3,5,7,9,11,13,15 -B 2 -z 4 -n -13 -d DBPSK -H 0.13,1.1 >> validate.out
echo "*** Message correlation BPSK" | tee -a validate.out
/opt/install/wsgc_test/bin/wsgc_test -s 4096 -c 1023 -C 10 -t 2.12 -r 0.0 -N 4 -I 0 -R 16 -P 1 -B 4 -F 31 -n -24 --cuda >> validate.out
echo "*** Training sequence correlation DBPSK" | tee -a validate.out
/opt/install/wsgc_test/bin/wsgc_test -s 4096 -c 1023 -C 140 -t 0.0 -r 0.0 -N 4 -I 0 -R 32 -Y 32 -M 64 -z 4 -d DBPSK -n -13 --simulate-trn >> validate.out
echo "*** Message correlation MFSK" | tee -a validate.out
#/opt/install/wsgc_test/bin/wsgc_test -s 4096 -d "MFSK:3,0,-32" -t 0.0 -M 64 -R 16 -n -20 -f "W:0.0,1.0,5.0,0.0" >> validate.out
/opt/install/wsgc_test/bin/wsgc_test -s 4096 -d MFSK:3,-1,-32 -t 0.0 -M 64 -R 63 -n -24 -H 1.75 >> validate.out
echo "*** MFSK with JT65 classic source encoding and Reed-Solomon with RSSoft library" | tee -a validate.out
/opt/install/wsgc_test/bin/wsgc_test -s 4096 -d "MFSK:3,0,-32" -t 0.0 -M 64 -R 16 -n -20 -f "W:0.0,1.0,5.0,0.0" -j "JT65::F6HTJ F4EXB JN33" -O RS/63,12,16,10,4:iarith,regex,.*F4EXB.* >> validate.out
echo "*** MFSK with Convolutional Coding with CCSoft library, Layland-Lushbaugh rate 1/2 k=32 code" | tee -a validate.out
/opt/install/wsgc_test/bin/wsgc_test -s 4096 -d "MFSK:3,0,-32" -t 0.0 -M 2 -R 36 -n -26 -O CC/32/4073739089,3831577671/fano:-1.5,-8,1,2000000,-4,20000000,-20 >> validate.out
echo "*** Message correlation WSGCE with CCSoft Convolutional Coding" | tee -a validate.out
/opt/install/wsgc_test/bin/wsgc_test -s 4096 -d WSGCE -t 0.0 -M 2 -R 16 -n -20 -O "CC/3/5,7/stack,,matchmsg:-1.4,6250000,-24,8,0.05" -B 4 >> validate.out
