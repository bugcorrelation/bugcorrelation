
% a axis --number of neighbours that are buggy
% X axis --all metrics
% y axis --probability that ego is buggy
datasets={'HTTPClient','Jackrabbit','Lucene','Rhino'};
graphtypes={'callgraph','hiegraph','cochgraph'};
reverse_ege_test=false;
results_LRC=zeros(4,24);
results_GEE=zeros(4,24);

for d=1:4
   coef=zeros(8,1);
   for g=[1:3]
      disp('=========================================================');  
      disp(['Data Set:',datasets{d},'  Graph Type:',graphtypes{g}]);
      disp(' ');
      filepath=['./data/',char(datasets{d}),'\'];
      bugstat_file=[filepath,'bugstat.txt'];
      bugstat0=load(bugstat_file);
      release_file=[filepath,'/bugs/release_date_revision.csv'];
      [releases,~,~]=textread(release_file,'%s%s%s','delimiter', ',');
      release=releases(1);
      graph_file=[filepath,char(graphtypes{1}),char(release),'.txt'];
      graph=load(graph_file);
      [filenum,~]=size(graph);
      versionnum=size(bugstat0,2)-1;
      %versionnum=6;
      bugstat=zeros(filenum,versionnum)-1;
      bugstat(bugstat0(:,1),:)=bugstat0(:,2:versionnum+1);
      metric_file=[filepath,'metrics',char(release),'.txt'];
      metrics0=load(metric_file);
      metricnum=size(metrics0,2)-1;

      metrics_bug_relation=zeros(filenum*(versionnum-1),metricnum);
      
      % for GEE
      nei_prev_bug=zeros(filenum*(versionnum-1),1);
      nei_curr_bug=zeros(filenum*(versionnum-1),1);
      ego_prev_bug=zeros(filenum*(versionnum-1),1);
      ego_curr_bug=zeros(filenum*(versionnum-1),1);
      id=zeros(filenum*(versionnum-1),1);
      t=zeros(filenum*(versionnum-1),1);

      for INFLUENCED_VERSION=2:versionnum
          
         %%graph construction
         release_curr=releases(INFLUENCED_VERSION);
         release_prev=releases(INFLUENCED_VERSION-1);

             graph_prev_file=[filepath,char(graphtypes{g}),char(release_prev),'.txt'];
             graph_prev=load(graph_prev_file);
             graph_curr_file=[filepath,char(graphtypes{g}),char(release_curr),'.txt'];
             graph_curr=load(graph_curr_file);
         
         
         if reverse_ege_test
             graph_prev=graph_prev';
             graph_curr=graph_curr';
         end
         
         
         metric_file=[filepath,'metrics',char(release_curr),'.txt'];
         metrics0=load(metric_file);
         metrics=zeros(filenum,size(metrics0,2)-1);
         metrics(metrics0(:,1),:)=metrics0(:,2:size(metrics0,2));
         
         metrics_bug_relation(filenum*(INFLUENCED_VERSION-2)+1:filenum*(INFLUENCED_VERSION-2)+filenum,:)=metrics(:,:);
         ego_prev_bug(filenum*(INFLUENCED_VERSION-2)+1:filenum*(INFLUENCED_VERSION-2)+filenum)=bugstat(:,INFLUENCED_VERSION-1);
         ego_curr_bug(filenum*(INFLUENCED_VERSION-2)+1:filenum*(INFLUENCED_VERSION-2)+filenum)=bugstat(:,INFLUENCED_VERSION);
         for i=1:filenum
            prev_nei_ind=(graph_prev(i,:)>0)';
            prev_exist_ind=bugstat(:,INFLUENCED_VERSION-1)>=0;
            prev_nei_exist_ind=prev_nei_ind & prev_exist_ind;
            prev_buggy_nei_ind=bugstat( prev_nei_exist_ind,INFLUENCED_VERSION-1)>0;
            nei_prev_bug(filenum*(INFLUENCED_VERSION-2)+i)=sum(prev_buggy_nei_ind);
            
            curr_nei_ind=(graph_curr(i,:)>0)';
            curr_exist_ind=bugstat(:,INFLUENCED_VERSION)>=0;
            curr_nei_exist_ind=curr_nei_ind & curr_exist_ind;
            curr_buggy_nei_ind=bugstat(curr_nei_exist_ind,INFLUENCED_VERSION)>0;
            nei_curr_bug((versionnum-1)*(i-1)+INFLUENCED_VERSION-1)=sum(curr_buggy_nei_ind);
            id(filenum*(INFLUENCED_VERSION-2)+i)=i;
            t(filenum*(INFLUENCED_VERSION-2)+i)=INFLUENCED_VERSION-1;
         end
      end
      
      nei_bug=nei_prev_bug+nei_curr_bug;
      
      
      %Logistic regression of buggy neighbor numbers
      disp('Logistic Regression for LRC...');
      for rr=2:100
         remain=sum(nei_prev_bug(:,1)>rr);
        av=(sum(nei_prev_bug))/(size(nei_prev_bug,1));
        if remain<av
           R=rr
           break;
        end
      end
      index=ego_curr_bug>=0 & nei_prev_bug<5;
      a=nei_prev_bug(index,1)+1;%neighbor previous bug number
      b=ego_prev_bug(index,1)+1;%ego previous bug number;
      c=nei_curr_bug(index,1)+1;%neighbor current bug number;
      noc=metrics_bug_relation(index,1)./2;%number of changes
      noa=metrics_bug_relation(index,2);%number of authors
     % nest=metrics_bug_relation(index,35);%number of max nesting
      loc=metrics_bug_relation(index,24)./20;%LOC is the 23th metric
      cyc=metrics_bug_relation(index,3)./5;% avg of cyclomatic
      ess=metrics_bug_relation(index,6)./2;% avg of essential complexity
       % X=[normalization([noc,noa,nest,cyc,ess,loc,b]),a];% merge the LOC and neighbor bug number as a new metric matrix

      M=[loc,noa,cyc,noc,c,b];
      X=[M,a]; 
      y=ego_curr_bug(index)>0;
      [alpha,~,stat]= glmfit(X,y,'binomial');
      nbn_index=size(X,2)+1;% the order index of the buggy neibor number in the result
      LRC=alpha(nbn_index);
      se=stat.se(nbn_index);
      p=stat.p(nbn_index);
      CI_l=alpha(nbn_index)-1.96*se;
      CI_h=alpha(nbn_index)+1.96*se;
      
      LRC=[LRC,CI_l,CI_h,p]
      results_LRC(d,((g-1)*4+1):((g-1)*4+4))=LRC;
      


      %%%%%%%%%%%%% GEE %%%%%%%%%%%%%
      disp('GEE...');
       R1=median(unique(nei_prev_bug))+1%
       R2=2;
      for rr=2:100
        remain=sum(nei_prev_bug(:,1)>rr);
        av=(sum(nei_prev_bug,1))/(R1);
        if remain<av
           R2=rr
           break;
        end
      end
      if g==9 index=(ego_curr_bug>=0) & (nei_prev_bug<3);
      else index=(ego_curr_bug>=0) & (nei_prev_bug<6);
      end
      
      noc=metrics_bug_relation(index,1);%number of changes
      noa=metrics_bug_relation(index,2);%number of authors
      loc=metrics_bug_relation(index,24);%LOC is the 23th metric
      cyc=metrics_bug_relation(index,37);% sum of cyclomatic
      ess=metrics_bug_relation(index,40);% sum of essential complexity
      tt=metrics_bug_relation(index,:);

      M=[noc,noa,cyc,ess,loc];
      y=ego_curr_bug(index)>0;
     % pwd;  currentFolder = pwd;  addpath([currentFolder,'/GEEQBOX/']);
      varnames={'number of changes','number of authors','Cyclomatic','Essent_Complex','LOC','prev_nei_bug','prev_ego_bug','constant'};
     % varnames=cell(1,size(M,2)+4);  varnames(size(M,2)+1:size(M,2)+4)={'neigh_curr_bug','nei_prev_bug','prev_ego_bug','constant'}; for dd=1:size(M,2)varnames(dd)={['metric',dd]};end
      X=[M,nei_prev_bug(index)+1,ego_prev_bug(index)+1,ones(sum(index),1)];
      [betahat, alphahat, results] = gee(id(index), y, t(index), X, 'b','equi',varnames);
      GEE0=results.model;
      GEE0=[GEE0{7,2},GEE0{7,6},GEE0{7,7},GEE0{7,5}];
      GEE=[results.model{8,2},results.model{8,6},results.model{8,7},results.model{8,5}];
      results_GEE(d,((g-1)*4+1):((g-1)*4+4))=GEE;

   end
end

csvwrite('./results_c1.csv',[results_LRC',results_GEE']');

