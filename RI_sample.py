#/usr/bin/python
# Sample source code for "Recurrent Infomax generates synfire chains,
# neuronal avalanches, and simple cell-like selectivity"
# (Neural Computation 21(4) 1038-1067) by Takuma Tanaka (tanaka.takuma@gmail.com).
# To reproduce Fig. 6
# $ python RI_sample.py -N 50 --pbar 0.05 --pmax 0.95 -B 2500 -S 0 -T 50000 --name sequence
#
# To reproduce Fig. 8
# $ python RI_sample.py -N 50 --pbar 0.01 --pmax 0.5 -B 20000 -S 5000 -T 10000 --name avalanche

import os;
import sys;
import math;
import numpy;
import random;
import argparse;

parser = argparse.ArgumentParser(description="Recurrent Infomax sample program");
parser.add_argument("-N", action="store", type=int, dest="N", default=20, help="number of neurons");
parser.add_argument("--pbar", action="store", type=float, dest="pbar", default=0.01, help="mean firing rate");
parser.add_argument("--pmax", action="store", type=float, dest="pmax", default=0.95, help="reliability of the neurons");
parser.add_argument("-B", action="store", type=int, dest="B", default=10000, help="number of learning blocks");
parser.add_argument("-S", action="store", type=int, dest="S", default=0, help="number of blocks counting the number of avalanches");
parser.add_argument("-T", action="store", type=int, dest="T", default=10000, help="steps in a block");
parser.add_argument("--epsilon", action="store", type=float, dest="epsilon", default=0.01, help="learning rate of h");
parser.add_argument("--eta", action="store", type=float, dest="eta", default=0.2, help="learning rate of W");
parser.add_argument("--initialskip", action="store", type=int, dest="initialskip", default=5, help="initial blocks without the update of W");
parser.add_argument("--loginterval", action="store", type=int, dest="loginterval", default=5, help="logging interval");
parser.add_argument("--name", action="store", type=str, dest="name", required=True, help="output file name");
parser.add_argument("--seed", action="store", type=int, dest="seed", default=1, help="seed of random number generators");

args = parser.parse_args();
N = args.N;
pbar = args.pbar;
pmax = args.pmax;
B = args.B;
S = args.S;
T = args.T;
loginterval = args.loginterval;
initialskip = args.initialskip;
epsilon = args.epsilon;
eta = args.eta;
name = args.name;
seed = args.seed;

random.seed(seed);
numpy.random.seed(seed);
meanI = 0;

x = numpy.array([0]*N);
W = numpy.random.rand(N,N)-0.5;
h = numpy.array([0.]*N);
p = numpy.array([pbar]*N);

histmax = 1000;
shist = [0]*histmax;

myfile = open("log_"+name+".txt", "w");

burstsize = 0;
raster = [0]*T;

def onestep(corrmode):
    global x, h, W, burstsize, shist, C, Eyx;
    y = (numpy.sign(pmax/(1+numpy.exp(-W.dot(x)+h))-numpy.random.rand(N))+1)/2;
    h += epsilon*(y-p);
    if corrmode:
        C += numpy.outer(y-p, y-p);
        Eyx += numpy.outer(y-p, x-p);
    x = y;
    if longiter>=B-S:
        if burstsize>0:
            if numpy.sum(x)==0:
                shist[burstsize] += 1;
                burstsize = 0;
            else:
                burstsize += int(numpy.sum(x));
        elif numpy.sum(x)>0:
            burstsize = int(numpy.sum(x));

for longiter in range(B):
    for shortiter in range(2*T):
        onestep(False);

    C = numpy.zeros((N,N));
    Eyx = numpy.zeros((N,N));
    for shortiter in range(T):
        onestep(True);
        if longiter==B-1:
            raster[shortiter] = x;
    C /= T;
    Eyx /= T;
    
    D = numpy.r_[numpy.c_[C, numpy.transpose(Eyx)], numpy.c_[Eyx, C]];
    C_ = numpy.linalg.inv(C);
    CC_ = numpy.r_[numpy.c_[C_, numpy.zeros((N, N))], numpy.c_[numpy.zeros((N, N)), C_]];
    D_ = numpy.linalg.inv(CC_.dot(D)).dot(CC_);
    if initialskip<=longiter and longiter<B-S:
        V = (numpy.ones((N, N))-numpy.eye(N))*(2*C_-D_[0:N,0:N]-D_[N:2*N,N:2*N]);
        D_t = numpy.transpose(D_);
        U = (numpy.outer(1-2*p, 1-2*p)*Eyx+numpy.outer(p*(1-p), p*(1-p))-Eyx*Eyx)*(D_[N:2*N,0:N]+D_t[N:2*N,0:N]);
        W += eta*(numpy.transpose(C).dot(V).dot(Eyx)+numpy.transpose(numpy.transpose(Eyx).dot(V).dot(C))-U)/2;
    meanI += -math.log(numpy.linalg.det(CC_.dot(D)))/(2*math.log(2.0));
    
    if longiter%loginterval==loginterval-1:
        message = str(longiter+1)+" "+str(meanI/loginterval)+"\n";
        sys.stdout.write(message)
        myfile.write(message)
        myfile.flush();
        meanI = 0;
myfile.close();

if S>0:
    myfile = open("hist_avalanche_"+name+".txt", "w");
    for i in range(histmax):
        if shist[i]>0:
            myfile.write(str(i)+" "+str(shist[i]/float(sum(shist)))+"\n");
    myfile.close();

braster = [0]*T;
myfile = open("raster_"+name+".txt", "w");
for i in range(T):
    bvec = "";
    for j in range(N):
        if raster[i][j]!=0:
            myfile.write(str(i)+" "+str(j)+"\n");
            bvec += "1";
        else:
            bvec += "0";
    braster[i] = bvec;
myfile.close();

repcount = [(x, braster.count(x)) for x in set(braster)];
rephist = [0]*N;
for x,num in repcount:
    rephist[x.count('1')] += num;
myfile = open("pattern_repetition_"+name+".txt", "w");
for i in range(1,N):
    if rephist[i]>0:
        myfile.write(str(i)+" "+str(rephist[i])+"\n");
myfile.close();

myfile = open("connection_"+name+".txt", "w");
for i in range(N):
    for j in range(N):
        myfile.write(str(W[i][j])+" ");
    myfile.write("\n");
myfile.close();
