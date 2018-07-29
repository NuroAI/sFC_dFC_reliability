
% sliding window
n_ROI=160;
mat=tril(ones(n_ROI,n_ROI),-1);
ind_mat=find(mat(:)>0);

for windowsize=10:10:200
    
    TR=0.72;
    % windowsize=40;  % in the unit of how many TR
    w=windowsize*TR;
    sigma=3;        % width of the gaussian kernel (in the unit of how many TR, used to convolve with rectanglar window to derive a tapered window)
    
    n=820;       % number of subjects
    nv=n_ROI*(n_ROI-1)/2;   % number of FC
    nr=4;       % number of runs
    
%     sw_d1=zeros(nv,n,nr);
%     sw_d2=sw_d1;
%     sw_d3=sw_d1;
    sw_d4=zeros(nv,n,nr);
    
    folder={'C:/Users/Chao/Downloads/820 common subjects/DOS160/run1/','C:/Users/Chao/Downloads/820 common subjects/DOS160/run2/','C:/Users/Chao/Downloads/820 common subjects/DOS160/run3/','C:/Users/Chao/Downloads/820 common subjects/DOS160/run4/'};
    
    for Run=1:4;
        
        path=folder{Run};
        flist=dir(path);
        flist={flist.name};
        flist=flist(3:n+2);
        
        %         f_loc=zeros(100,12720);
        
        for i=1:n
            
            load(strcat(path,flist{i}))
            tc=tc(:,dos_new);
            temp=tc;
            %     temp=tapered_sliding_window(tc, windowsize, sigma);
            temp=y_IdealFilter(temp, 0.72, [0.01 0.1]);
            temp=y_IdealFilter(temp, 0.72, [1/w 0]);  % high-pass to avoid spurious fluctions
            
            temp=sliding_window(temp, windowsize);
            
            temp=reshape(temp,n_ROI*n_ROI,1200);
            temp=temp(ind_mat,:);
            temp(:,isnan(temp(1,:)))=[];
            
%             sw_d1(:,i,Run)=mean(temp,2);
%             sw_d2(:,i,Run)=std(temp,[],2);
%             [sw_d3(:,i,Run)]=alff_2D(temp,w);
            
            tc_corr=temp';
            
            temp1=tc_corr-repmat(median(tc_corr),size(tc_corr,1),1);
            temp2=sign(temp1);
            for j=1:nv
                temp3=abs(temp1(:,j));
                temp4=temp2(:,j);
                ind=[1; find((temp4(2:end)~=temp4(1:end-1)).*(temp4(2:end)~=0))+1];  % cross median point locations starting point
                l=length(ind);
                %excursion_startend=[ind(1:end-1)-double(temp4(ind(1:end-1)-1)==0) ind(2:end)-1];
                excursion_startend=[ind [ind(2:end)-1;size(tc_corr,1)]];
                excursion_len=(excursion_startend(:,2)-excursion_startend(:,1)+1);
                excursion_hei=zeros(l,1);
                for k=1:l, excursion_hei(k)=max(temp3(excursion_startend(k,1):excursion_startend(k,2))); end
                sw_d4(j,i,Run)=sum((excursion_len.^(0.9)).*excursion_hei);
            end
                        
            %
            %             temp=alff_2D_location(temp,w);
            %             f_loc(i,:)=temp;
            %             pl=zeros(160*160,1);
            %             % pl(ind_mat)=mean(Ct2,2);
            %             % pl(ind_mat)=std(Ct2,[],2);
            %             pl(ind_mat)=temp;
            %             pl=reshape(pl,160,160);
            %             pl=pl+pl';
            %             figure;imagesc(pl)
            %             axis square
            %             colorbar
            i
        end
        
        %         for j=1:10
        %         figure;hist(f_loc(:,j))
        %         end
    end
    
%     sw_d1=permute(sw_d1,[3,2,1]);
%     sw_d2=permute(sw_d2,[3,2,1]);
%     sw_d3=permute(sw_d3,[3,2,1]);
    sw_d4=permute(sw_d4,[3,2,1]);
    
    sn=strcat('dynamic_820_slidingwindow_bandpass_highpass_win',num2str(windowsize,'%03d') ,'_leonardi2015.mat');
%     save(sn,'sw_d1','sw_d2','sw_d3','sw_d4');
    save(sn,'sw_d4','-append')
    %     save(sn,'sw_d2');
    windowsize
end




%% static FC
TR=0.72;
windowsize=40;  % in the unit of how many TR
w=windowsize*TR;
sigma=3;        % width of the gaussian kernel (in the unit of how many TR, used to convolve with rectanglar window to derive a tapered window)

n_ROI=160;
mat=tril(ones(n_ROI,n_ROI),-1);
ind_mat=find(mat(:)>0);
nv=length(ind_mat);

n=820;       % number of subjects
nr=4;       % number of runs

%     sw_d1=zeros(nv,n,nr);
%     sw_d2=sw_d1;
%     sw_d3=sw_d1;
sw_d4=zeros(nv,n,nr);

X0=zeros(nr,n,nv);

folder={'C:/Users/Chao/Downloads/820 common subjects/DOS160/run1/','C:/Users/Chao/Downloads/820 common subjects/DOS160/run2/','C:/Users/Chao/Downloads/820 common subjects/DOS160/run3/','C:/Users/Chao/Downloads/820 common subjects/DOS160/run4/'};
load('E:\imp_Chao\projects\07_TRT\dos_neworder.mat')

for Run=1:4;
    
    path=folder{Run};
    flist=dir(path);
    flist={flist.name};
    flist=flist(3:n+2);
    
    for i=1:n
        
        load(strcat(path,flist{i}))
        tc=tc(:,dos_new);
        temp=tc;
        %     temp=tapered_sliding_window(tc, windowsize, sigma);
        temp=y_IdealFilter(temp, TR, [0.01 0.1]);
        temp=corr(temp);
        X0(Run,i,:)=temp(ind_mat);
        i
    end
        
        
end

        
