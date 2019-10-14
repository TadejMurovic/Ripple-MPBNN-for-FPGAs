clear all; close all;
LT    = 192;
SEL_T = 1;

for la = 0:1
    
    % Config
    if la == 0
        N = LT; % Neurons per hidden layer
        M = 593;
        sel = SEL_T; % sel = 1 Ripple Arch; sel = 0 Straight-Forward MPBNN code generation
    end
    if la == 1
        N = 1; % Neurons per hidden layer
        M = LT;
        sel = 0; % sel = 1 Ripple Arch; sel = 0 Straight-Forward MPBNN code generation
    end
    pathw = strcat('Net_Params\dump_w_gen',num2str(LT),'_',num2str(la),'.txt');
    patht = strcat('Net_Params\dump_t_gen',num2str(LT),'_',num2str(la),'.txt');



    %--------------------------------------------------------------------------
    % Read Data
    fileID = fopen(pathw,'r');
    formatSpec = '%d';
    cols = M;
    rows = N;
    w = fscanf(fileID,formatSpec,[cols rows])';
    fclose(fileID);
    fileID = fopen(patht,'r');
    formatSpec = '%d';
    th = fscanf(fileID,formatSpec,[rows 1])';
    fclose(fileID);   

    if sel == 1
        % Optimize Order
        ordered = hamming_genetic(w,5000,1000);
        word = w(ordered,:);
    end
 
    if sel == 1
        n = N;
        m = M;
        wts = word';
        a0 = [];
        a1 = [];
        d0 = [];
        d1 = [];
        for ni = n 
            counter = 0;
            dep = 0;
            for i = 0:log2(ni)-1
               dep = dep + 1;
               for j = 0:ni-1 
                  if j >= 2^i
                    counter = counter + 1;    
                  else
                    counter = counter + 0;   
                  end
               end
            end
            a0 = [a0 counter];
            d0 = [d0 dep];

            counter = 0;
            dep = 0;
            for d = 0:log2(ni)-1
                dep = dep + 1;
                for i = 0:(ni/2^d - 2) 
                    counter = counter + 1;       
                end
            end  
            for d = (log2(ni)-1):-1:0
                dep = dep + 1;
                for i = 0:(ni/2^d - 1)  
                    if i > 0
                       if (mod(i,2) ~= 0)
                           counter = counter + 0;
                       else
                           counter = counter + 1;
                       end
                    else
                       counter = counter + 0; 
                    end
                end
            end
            a1 = [a1 counter];
            d1 = [d1 dep];  
        end
        levels = dep;

        cat_path = strcat('model_',num2str(la),'.txt'); 
        fid = fopen(cat_path,'w');

        % Generate
        fprintf(fid,'`timescale 1ns / 1ps');
        fprintf(fid,'\n\n');
        fprintf(fid,strcat('module layer_',num2str(la),'(in, out'));
        fprintf(fid,');\n\n');

        fprintf(fid,'  input [%d:0] in;\n',m-1);
        fprintf(fid,'  output reg [%d:0] out;\n',n-1);
        fprintf(fid,'\n');
            bits = ceil(log2(m)+1); 
            bts = ceil(log2(m*n)+1);
            bts = bts + 1;
            bits = bts;
            %levels = ceil(log2(n+1));
        fprintf(fid,'  wire [%d:0] signal [%d:0];\n',m-1,n-1);
        fprintf(fid,'  wire [%d:0] difsum [%d:0];\n',bits-1,n-1);
        fprintf(fid,'  wire [%d:0] cumsum [%d:0][%d:0];\n',bts-1,levels-1,n-1);
        fprintf(fid,'  wire [%d:0] setsum;\n',bits-1);
        fprintf(fid,'  wire signed [%d:0] eval [%d:0];\n',bts,n-1);
        fprintf(fid,'\n\n');
        fprintf(fid,'  assign signal[0] = in;\n');
        fprintf(fid,'  assign setsum    = ');
        for j = 1:m
            fprintf(fid,'signal[0][%d] ',j-1);
            if j == m
                fprintf(fid,';\n\n'); 
            else
                fprintf(fid,'+ ');
            end
        end
        wts = [ones(m,1) wts];

        for i = 1:n
            mask = xor(wts(:,i+1),wts(:,i));
            mask = mask';
            mask = fliplr(mask);
            s(i) = sum(mask); 
            if i < n
                fprintf(fid,'  assign signal[%d] = signal[%d] ^ %d''b',i,i-1,m); 
                for j = 1:m
                    fprintf(fid,'%d',mask(j));    
                end
            end
            mask = fliplr(mask);
            fprintf(fid,';\n');
            fprintf(fid,'  assign difsum[%d] = ',i-1); 
            for j = 1:m
                chk = 0;
                if mask(j) == 1
                    fprintf(fid,'signal[%d][%d] ',i-1,j-1);
                    chk = 1;
                end
                if j == m || sum(mask(j+1:end)) == 0 
                    if j == 1
                      fprintf(fid,'1''b0');     
                    end
                    fprintf(fid,';\n'); 
                    break;
                else
                    if chk == 1
                        fprintf(fid,'+ ');
                    end
                end
            end
        end
        fprintf(fid,'\n\n');

        for i = 0:n-1
            fprintf(fid,'  assign cumsum[%d][%d] = difsum[%d];\n',0,i,i);  
        end
        fprintf(fid,'\n\n');

        counter = 0;
        dep = 0;
        ni = n;
        for d = 0:log2(ni)-1
            dep = dep + 1;
            seti = 2^(d+1)-1:2^(d+1):ni;
            for i = 0:ni-1
                if sum(i==seti)==1
                    fprintf(fid,'  assign cumsum[%d][%d] = cumsum[%d][%d] + cumsum[%d][%d];\n',dep,i,dep-1,i,dep-1,i-2^(d));  
                    counter = counter + 1; 
                else
                    fprintf(fid,'  assign cumsum[%d][%d] = cumsum[%d][%d];\n',dep,i,dep-1,i);  
                end
            end
        end  
        for d = ceil(log2(ni))-2:-1:0
            dep = dep + 1;
            seti = 2^(d+1):2^(d+1):ni;
            seti(seti > ni-1) = [];
            seti = seti + 2^d-1;
            for i = 0:ni-1
                if sum(i==seti)==1
                   fprintf(fid,'  assign cumsum[%d][%d] = cumsum[%d][%d] + cumsum[%d][%d];\n',dep,i,dep-1,i,dep-1,i-2^(d));  
                   counter = counter + 1; 
                else
                   fprintf(fid,'  assign cumsum[%d][%d] = cumsum[%d][%d];\n',dep,i,dep-1,i);  
                end
            end
        end
        fprintf(fid,'\n\n');

        for i = 0:n-1
            fprintf(fid,'  assign eval[%d] = $signed(cumsum[%d][%d] << 1''b1 );\n',i,levels-1,i);  
        end
        fprintf(fid,'\n\n');

        os = s;
        oth = th;
        s = cumsum(s);

        fprintf(fid,'  always @* begin\n');
        for j = 1:n
            thresh = th(j)-s(j);    
            bit    = bts;
            if thresh >= 0  
                fprintf(fid,'    out[%d]= (eval[%d]<=($signed(setsum)-$signed(%d''sd%d)));\n',j-1,j-1,bit+1,abs(thresh));
            else
                fprintf(fid,'    out[%d]= (eval[%d]<=($signed(setsum)-$signed(-%d''sd%d)));\n',j-1,j-1,bit+1,abs(thresh));
            end
        end
        fprintf(fid,'  end\n');
        fprintf(fid,'\n\n'); 

        fprintf(fid,'\n\n');  
        fprintf(fid,'endmodule');
        fclose(fid);
    else  
      layer_s = LT;
      for i = 1:1
        % Make standard model verilog code
        Nn = M;
        M = N;
        N = Nn;
        
        lay = i;
        bits = ceil(log2(N+1)); 
        cat_path = strcat('model_',num2str(la),'.txt'); 
        fid = fopen(cat_path,'w');

        % Generate
        fprintf(fid,'`timescale 1ns / 1ps');
        fprintf(fid,'\n\n');
        fprintf(fid,strcat('module layer_',num2str(la),'(in, out'));
        fprintf(fid,');\n\n');

        fprintf(fid,'  input [%d:0] in;\n',N-1);
        fprintf(fid,'  output reg [%d:0] out;\n',M-1);
        for ii = 0:M-1
           %fprintf(fid,'  output reg out%d;\n',i);
           fprintf(fid,'  reg [%d:0] t%d;\n',bits-1,ii);
           fprintf(fid,strcat('  reg [%d:0] w%d = %d''b'),N-1,ii,N);
           for j = 1:N
               fprintf(fid,'%d',w(ii+1,j)); 
           end
           fprintf(fid,';\n');
           fprintf(fid,'  reg [%d:0] th%d = %d''d%d;\n',bits-1,ii,bits,abs(th(ii+1)));
           fprintf(fid,'  reg [%d:0] weighted%d;\n',N-1,ii);
        end

        fprintf(fid,'\n\n');
        fprintf(fid,'  integer idx;\n\n');
        fprintf(fid, '  always @* begin\n');
        fprintf(fid, '    for( idx = 0; idx<%d; idx = idx + 1) begin\n',N);
        for ii = 0:M-1
           fprintf(fid,'      weighted%d[idx] = ((w%d[idx])~^(in[idx]));\n',ii,ii); 
        end
        fprintf(fid,'    end\n');
        fprintf(fid,'  end\n');

        fprintf(fid,'\n\n');
        fprintf(fid, '  always @* begin\n');
        for ii = 0:M-1
           fprintf(fid,'    t%d = ',ii); 
           for j = 0:N-2
              fprintf(fid, 'weighted%d[%d] + ',ii,j); 
           end
           j = j + 1;
           fprintf(fid, 'weighted%d[%d];\n',ii,j);
        end
        fprintf(fid,'  end \n');

        fprintf(fid,'\n\n');
        fprintf(fid, '  always @* begin\n');
        for ii = 0:M-1
           fprintf(fid,'    out[%d] = t%d >= th%d;\n',ii,ii,ii); 
        end    
        fprintf(fid,'  end \n');

        fprintf(fid,'endmodule');
        fclose(fid);  
      end
    end
end
fclose('all')
return
