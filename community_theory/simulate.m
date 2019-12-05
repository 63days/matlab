function test()
bit_size = 1000;
symbol_size = bit_size/2;
bArr = binornd(1, 1/2, bit_size, 1); % source coding complete.
sArr = zeros(symbol_size, 1);
u = 1/sqrt(2);
s1 = complex(u, u);
s2 = complex(-u, u);
s3 = complex(-u, -u);
s4 = complex(u, -u);
for i = 1 : 2 : bit_size-1 %constellation
    if(bArr(i) == 0 && bArr(i + 1) == 0) %(00)
        sArr((i+1)/2) = s1;
    elseif(bArr(i) == 0 && bArr(i + 1) == 1) %(01)
        sArr((i+1)/2) = s2;
    elseif(bArr(i) == 1 && bArr(i + 1) == 0) %(10)
        sArr((i+1)/2) = s3;
    else                                   %(11)
        sArr((i+1)/2) = s4;
    end
end
%%%%%%%%%%%%%%%Constellation mapping%%%%%%%%%%%%%%
subplot(2,2,1);
plot(sArr, 'o');
title('Constellation');
xlabel('Re');
ylabel('Im');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
T = 1/15000; %bit rate = 30kbps 가정
t = -symbol_size*T : 0.01*T : symbol_size*T;
ts = downsample(t, 100);

x = zeros(1, numel(t));

for i = 1 : symbol_size   %x(t)function generation
    x = x + sArr(i)*gtx(t-i*T);
end
xs = downsample(x, 100);

%%%%%%%%%%%%%%x(t) t-domain plot%%%%%%%%%%%%%
subplot(2, 2, 2);
plot(t, x); hold on;
plot(ts, xs, 'o');
title('x(t) time-domain');
xlabel('t');
ylabel('x(t)');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

v = normrnd(0, 0.001, [1 length(t)]); %AWGN channel

fc = 2000*10^6; %carrier freq
xp = real(x.*exp(1i*2*pi*fc*t));

%%%%%%%%%%%%%X(f), Xp(f) psd%%%%%%%%%%%%%%%%%
subplot(2,2,3);
[p1 f1] = pspectrum(x);
plot(f1, pow2db(p1));
[p2 f2] = pspectrum(xp);
hold on;
plot(f2, pow2db(p2));
title('PSD');
xlabel('f');
ylabel('dB');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
yp = xp + v;
c = lowpass(yp.*cos(2*pi*fc.*t), 1/T, 100/T);
s = lowpass(yp.*sin(2*pi*fc.*t), 1/T, 100/T);
y = 2*c -1i*2*s; %baseband equivalent y(t)
yg = conv(y, gtx(t), 'same'); %y(t) after g_rx function
yg = yg/98.5540; %normalizing

ys = downsample(yg, 100);
ys = ys(symbol_size+2:2*symbol_size+1);
yd = zeros(1, symbol_size); %detecting y
for i = 1 : length(ys)
    temp = [abs(ys(i)-s1) abs(ys(i)-s2) abs(ys(i)-s3) abs(ys(i)-s4)];
    [minV index] = min(temp);
    if(index == 1)
        yd(i) = s1;
    elseif(index == 2)
        yd(i) = s2;
    elseif(index == 3)
        yd(i) = s3;
    else
        yd(i) = s4;
    end
end
k = 0;
for i = 1:symbol_size
    if(yd(i) == sArr(i))
        fprintf("%d 에서 같음 \n", i);
        k = k + 1;
    end
end
fprintf("number of successfully detected: %d \n", k);
%%%%%%%%%%%%%y mapping%%%%%%%%%%%%%%%%%
subplot(2, 2, 4);
plot(ys, 'o'); hold on;
plot(yd, '*');
title('y mapping');
xlabel('Re');
ylabel('Im');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

function y = gtx(t)
    T = 1/15000;
    y = zeros(1, numel(t));
    for i = 1:numel(t)
    if(abs(t(i)) < 7*T)
        y(i) = sinc(t(i)/T);
    else
        y(i) = 0;
    end
    end
end