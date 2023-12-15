function RI_sample(mode)

% Sample source code for "Recurrent Infomax generates synfire chains,
% neuronal avalanches, and simple cell-like selectivity"
% (Neural Computation 21(4) 1038-1067) by Takuma Tanaka.
%
% This source code runs on Matlab and Octave.
% If you find any bugs in this program, or have any question
% I would appreciate the tip
% (ttakuma@mbs.med.kyoto-u.ac.jp).
%
% To run this program:
%
% 1. If you are using Matlab, delete the line of code 'fflush(stdout);' 
% appearing later on in this file.
%
% 2. Go into the directory where the file RI_sample.m is located, and enter
%      matlab
% or
%      octave
% in the console.
%
% 3. Enter
%      RI_sample(1)
% for avalanche before learning (Fig. 8, N=50, p_max=0.5, \bar{p}_i=0.01),
%      RI_sample(2)
% for avalanche after learning (Fig. 8, N=50, p_max=0.5, \bar{p}_i=0.01),
%      RI_sample(3)
% for sequences before learning (Fig. 6, N=50, p_max=0.95, \bar{p}_i=0.05), or
%      RI_sample(4)
% for sequences after learning (Fig. 6, N=50, p_max=0.95, \bar{p}_i=0.05),
% and then the simulation will run.
%
% The simulation takes several hours to days in modes 2 and 4.
%
% The block number and approximate mutual information are displayed in
% the console and also output to the log files named 
% 'log_avalanche_before.txt' etc.
% At the end of the simulation, the histogram of the burst size is output to
% the files 'hist_avalanche_before.txt' and 'hist_avalanche_after.txt'
% in modes 1 and 2, respectively.
% Firing states of the neurons in the last block are output to the files
% 'raster_sequence_before.txt' and 'raster_sequence_after.txt'
% in modes 3 and 4, respectively.
% The numbers of repeated patterns are output to the files
% 'number_sequence_before.txt' and 'number_sequence_after.txt'
% in modes 3 and 4, respectively.
% The sequence size is output in the first column, and the number of
% occurrences is output in the second column.
% Optimized connection weights are output to the files named
% 'connection_avalanche_before.txt' etc.

if (mode~=1 & mode~=2 & mode~=3 & mode~=4)
  fprintf(1, 'Invalid mode.\n  Mode 1: avalanche before learning\n  Mode 2: avalanche after learning\n  Mode 3: sequence before learning\n  Mode 4: sequence after learning\n');
  return;
end;

rand('seed', 1);
meanI = 0;
logstep = 5;
factor = 100;

N = 50;
x = zeros(N, 1);
W = rand(N,N)-0.5;
h = zeros(N, 1);
epsilon = 0.01;
eta = 0.2;

histmax = 1000;
shist = zeros(histmax, 1);

if (mode==1 | mode==2)
  p = 0.01*ones(N, 1);
  pmax = 0.5;
  T = 10000;
  if (mode==1)
    fname = 'log_avalanche_before.txt';
    STAT = 5;
    BLOCKS = 5005;
  else
    fname = 'log_avalanche_after.txt';
    STAT = 15000;
    BLOCKS = 20000;
  end
else
  STAT = 5000;
  p = 0.05*ones(N, 1);
  pmax = 0.95;
  T = 50000;
  if (mode==3)
    fname = 'log_sequence_before.txt';
    BLOCKS = 5;
  else
    fname = 'log_sequence_after.txt';
    BLOCKS = 2500;
  end;
  raster = zeros(T, N);
end

myfile = fopen(fname, 'w', 'ieee-le');

state = 0;

for longiter=1:BLOCKS
  for shortiter=1:T
    y = pmax ./ (1 + exp(-W*x+h));
    y = rand(N, 1)<y;
    h = h+epsilon*(y-p);
    x = y;
    
    if (longiter>=STAT)
      if (state>0)
        if (sum(x)==0)
          shist(state) = shist(state)+1;
          state = 0;
        else
          state = state+sum(x);
        end
      elseif (sum(x)>0)
        state = sum(x);
      end;
    end;
  end;

  C = zeros(N,N);
  Eyx = zeros(N,N);

  for shortiter=1:T
    y = pmax ./ (1 + exp(-W*x+h));
    y = rand(N, 1)<y;
    h = h+epsilon*(y-p);
    C = C+(y-p)*(y-p)';
    Eyx = Eyx+(y-p)*(x-p)';
    x = y;
    if (longiter>=STAT)
      if (state>0)
        if (sum(x)==0)
          shist(state) = shist(state)+1;
          state = 0;
        else
          state = state+sum(x);
        end
      elseif (sum(x)>0)
        state = sum(x);
      end;
    elseif (longiter==BLOCKS & mode>=3)
      raster(shortiter,:) = x';
    end;
  end;

  C = C/T;
  Eyx = Eyx/T;
  D = [C Eyx' ; Eyx C];
  
  if (5<longiter & longiter<STAT)
    C_ = (factor*C)'^-1;
    CC_ = [C_ zeros(N) ; zeros(N) C_];
    D_ = ((factor*D*CC_)'^(-1))*CC_;
    V = (ones(N, N)-eye(N)) .* (2*C_-D_(1:N,1:N)-D_(N+1:2*N,N+1:2*N))*factor;
    D_t = D_';
    U = (((1-2*p)*(1-2*p)').*Eyx+(p.*(1-p))*(p.*(1-p))'-Eyx.*Eyx).*((D_(N+1:2*N,1:N))+(D_t(N+1:2*N,1:N)))*factor;
    W = W+eta*(C'*V*Eyx+(Eyx'*V*C)'-U)/2;
  end;

  I = (log(det(factor*C))-log(det(factor*D))/2)/log(2);
  meanI = meanI+I;
  
  fprintf(1, '%d ', longiter);
  disp(meanI/(mod(longiter-1, logstep)+1));

  fflush(stdout);  % Delete this line if you are using Matlab
  
  if (mod(longiter, logstep)==0)
    fprintf(myfile, '%d %e\n', longiter, meanI/logstep);
    if (mod(longiter, 100)==0)
      fclose(myfile);
      myfile = fopen(fname, 'a', 'ieee-le');
    end;
    meanI = 0;
  end;
  
end;

fclose(myfile);

if (mode<=2)
  if (mode==1)
    myfile2 = fopen('hist_avalanche_before.txt', 'w', 'ieee-le');
  else
    myfile2 = fopen('hist_avalanche_after.txt', 'w', 'ieee-le');
  end
  for i=1:histmax
    if (shist(i)>0)
      fprintf(myfile2, '%d %e\n', i, shist(i)/sum(shist));
    end;
  end;
else
  raster = [raster zeros(T, 2)];
  if (mode==3)
    myfile2 = fopen('raster_sequence_before.txt', 'w', 'ieee-le');
  else
    myfile2 = fopen('raster_sequence_after.txt', 'w', 'ieee-le');
  end;
  for i=1:T
    index = 0;
    for j=1:N
      fprintf(myfile2, '%d ', raster(i,j));
      index = index*2+raster(i,j);
    end;
    fprintf(myfile2, '\n');
    raster(i, N+1) = sum(raster(i, :));
    raster(i, N+2) = index;
  end;
  fclose(myfile2);
  [s, perm] = sort(raster(:,N+2));
  raster = raster(perm, :);
  [s, perm] = sort(raster(:,N+1));
  raster = raster(perm, :);

  repeated = zeros(N, 1);
  if (raster(1,N+2)==raster(2,N+2) & raster(2,N+1)>0)
    repeated(raster(1, N+1), 1) = repeated(raster(1, N+1), 1)+1;
  end;
  for i=2:T-1
    if (raster(i,N+1)>0 & (raster(i,N+2)==raster(i-1,N+2) | raster(i,N+2)==raster(i+1,N+2)))
      repeated(raster(i, N+1), 1) = repeated(raster(i, N+1), 1)+1;
    end;
  end;
  if (raster(T,N+1)>0 & raster(T,N+2)==raster(T-1,N+2))
    repeated(raster(T, N+1), 1) = repeated(raster(T, N+1), 1)+1;
  end;

  if (mode==3)
    myfile2 = fopen('number_sequence_before.txt', 'w', 'ieee-le');
  else
    myfile2 = fopen('number_sequence_after.txt', 'w', 'ieee-le');
  end;
  for i=1:N
    if (repeated(i)>0)
      fprintf(myfile2, '%d %d\n', i, repeated(i));
    end;
  end;
end;

fclose(myfile2);

if (mode==1)
  myfile3 = fopen('connection_avalanche_before.txt', 'w', 'ieee-le');
elseif (mode==2)
  myfile3 = fopen('connection_avalanche_after.txt', 'w', 'ieee-le');
elseif (mode==3)
  myfile3 = fopen('connection_sequence_before.txt', 'w', 'ieee-le');
else
  myfile3 = fopen('connection_sequence_after.txt', 'w', 'ieee-le');
end;

for i=1:N
  for j=1:N
    fprintf(myfile3, '%e ', W(i,j));
  end;
  fprintf(myfile3, '\n');
end;

fclose(myfile3);

