% a axis --number of neighbors that are buggy
% X axis --metrics
% y axis --probability that ego is buggy
clear; 

datasets={'HTTPClient','Jackrabbit','Lucene','Rhino'};
graphtypes={'callgraph','hiegraph','cochgraph','hybridgraph'};

num_train_releases=[2,6,2,3];% number of releases for training
num_basic_graph=3;% we have three basic graphs, call graph, hierarchy graph, and co-change graph.

for d=1:4
     disp('=========================================================');  
     disp(['Data Set:',datasets{d}]);
     disp(' ');
     filepath=['./data/',char(datasets{d}),'\'];
     bug_stat0=load([filepath,'bugstat.txt']);
     [releases,~,~]=textread([filepath,'/bugs/release_date_revision.csv'],'%s%s%s','delimiter', ',');
     graph_file=[filepath,char(graphtypes{1}),char(releases(1)),'.txt'];
     graph=load(graph_file);
     filenum=size(graph,1);
     release_num=size(releases,1)-1;
    
     metrics0=load([filepath,'metrics',char(releases(1)),'.txt']);
     metricnum=size(metrics0,2)-1;
     
     bug_stat=zeros(filenum,release_num)-1;
     bug_stat(bug_stat0(:,1),:)=bug_stat0(:,2:release_num+1);
     
     % for GEE
     metrics_stat=zeros(filenum*(release_num-1),metricnum);
     nei_prev_bug_basic=zeros(filenum*(release_num-1),8);% watched variable
     nei_prev_bug_hybrid=zeros(filenum*(release_num-1),1); % watched variable
     %nei_curr_bug=zeros(filenum*(release_num-1),num_basic_graph);
     ego_prev_bug=zeros(filenum*(release_num-1),1);
     ego_curr_bug=zeros(filenum*(release_num-1),1);
     id=zeros(filenum*(release_num-1),1);
     t=zeros(filenum*(release_num-1),1);
      
         
     disp('Loading data...');
     for INFLUENCED_VERSION=2:release_num
        INFLUENCED_VERSION
        release_curr=releases(INFLUENCED_VERSION);
        release_prev=releases(INFLUENCED_VERSION-1);
         
        metric_file=[filepath,'metrics',char(release_curr),'.txt'];
        metrics0=load(metric_file);
        metrics=zeros(filenum,size(metrics0,2)-1);
        metrics(metrics0(:,1),:)=metrics0(:,2:size(metrics0,2));
         
        metrics_stat(filenum*(INFLUENCED_VERSION-2)+1:filenum*(INFLUENCED_VERSION-2)+filenum,:)=metrics(:,:);
        ego_prev_bug(filenum*(INFLUENCED_VERSION-2)+1:filenum*(INFLUENCED_VERSION-2)+filenum)=bug_stat(:,INFLUENCED_VERSION-1);
        ego_curr_bug(filenum*(INFLUENCED_VERSION-2)+1:filenum*(INFLUENCED_VERSION-2)+filenum)=bug_stat(:,INFLUENCED_VERSION);
        id(filenum*(INFLUENCED_VERSION-2)+1:filenum*(INFLUENCED_VERSION-2)+filenum)=1:filenum;
        t(filenum*(INFLUENCED_VERSION-2)+1:filenum*(INFLUENCED_VERSION-2)+filenum)=INFLUENCED_VERSION-1;
         
        for g=[1:3]
         graph_prev_file=[filepath,char(graphtypes{g}),char(release_prev),'.txt'];
         graph_prev=load(graph_prev_file);
  

         for i=1:filenum
            prev_nei_ind=(graph_prev(i,:)>0)';
            prev_exist_ind=bug_stat(:,INFLUENCED_VERSION-1)>=0;
            prev_nei_exist_ind=prev_nei_ind & prev_exist_ind;
            prev_buggy_nei_ind=bug_stat( prev_nei_exist_ind,INFLUENCED_VERSION-1)>0;
            nei_prev_bug_basic((filenum*(INFLUENCED_VERSION-2)+i),g)=sum(prev_buggy_nei_ind);           
         end
        end
     end
    
      noc=metrics_stat(:,1);%number of changes
      noa=metrics_stat(:,2);%number of authors
      loc=metrics_stat(:,24);%LOC is the 23th metric
      cyc=metrics_stat(:,37);% sum of cyclomatic
      ess=metrics_stat(:,40);% sum of essential complexity
      M=[noc,noa,cyc,ess,loc];
      y=ego_curr_bug>0;
      
      
      
      
      
      
      
     %%%%%%%%%%%%%%%% Train LRCs %%%%%%%%%%%%%%%%%%%%
     disp('Training LRCs');
     lrc=zeros(3,1)+0.00001;
     for g=[1:num_basic_graph]
         disp(['Training...',datasets{d},'  Graph:',graphtypes{g}]);
         X=[M,nei_prev_bug_basic(:,g)+1, ego_prev_bug+1,ones(filenum*(release_num-1),1)];
         %%%%%%%%%%%%% GEE %%%%%%%%%%%%%
         disp('GEE...');
         index_train=zeros(filenum*(release_num-1),1);
         index_train(1:filenum*num_train_releases(d),:)=1;
         index_train=index_train>0;
         index= index_train &  (ego_curr_bug>=0) & (nei_prev_bug_basic(:,g)<6);
         varnames={'number of changes','number of authors','Cyclomatic','Essent_Complex','LOC','prev_nei_bug','prev_ego_bug','constant'};
         [betahat, alphahat, results] = gee(id(index), y(index), t(index), X(index,:), 'b','equi',varnames);
         GEE=[results.model{8,2},results.model{8,6},results.model{8,7},results.model{8,5}];
         lrc(g)=GEE(1);
     end
     
     for g=1:3
         if lrc(g)+max(lrc)<0 
             lrc(g)=0.00001;
         end
     end
     
      
 

      
      
      
      
       %%%%%%%%%%%%%%%% Test HybridGraph %%%%%%%%%%%%%%%%%%%
    for g=[1:4]
       disp(['Testing...',datasets{d},'  Graph:',graphtypes{g}]);
       if(g==4)% hybrid graph
          for INFLUENCED_VERSION=2:release_num
              INFLUENCED_VERSION
             release_curr=releases(INFLUENCED_VERSION);
             release_prev=releases(INFLUENCED_VERSION-1);
             graph_prev=zeros(filenum,filenum);
             graph_curr=zeros(filenum,filenum);
             for gg=[1:3]
                 graph_prev_file=[filepath,char(graphtypes{gg}),char(release_prev),'.txt'];
                 graph=load(graph_prev_file);
                 graph=graph>0;
                 if gg==2
                     graph_hie=graph;
                     %sum(graph_hie)
                 end
                 graph_prev=graph_prev+lrc(gg)*graph;
 
             end
             graph_prev=(graph_prev-max(lrc))>0;
            % sum(graph_prev)
             for i=1:filenum
                prev_nei_ind=(graph_prev(i,:)>0)';
                prev_exist_ind=bug_stat(:,INFLUENCED_VERSION-1)>=0;
                prev_nei_exist_ind=prev_nei_ind & prev_exist_ind;
                prev_buggy_nei_ind=bug_stat( prev_nei_exist_ind,INFLUENCED_VERSION-1)>0;
                nei_prev_bug_hybrid(filenum*(INFLUENCED_VERSION-2)+i)=sum(prev_buggy_nei_ind);
             end
          end
       end
    

      %%%%%%%%%%%%% GEE %%%%%%%%%%%%%
      disp('GEE...');
      index_test=ones(filenum*(release_num-1),1);
      index_test(1:filenum*num_train_releases(d),:)=0;
      index_test=index_test>0;
      if g==9
          X=[M,nei_prev_bug_hybrid, ego_prev_bug,ones(filenum*(release_num-1),1)];
          index=index_test & (ego_curr_bug>=0) & (nei_prev_bug_hybrid<3);
          
      else
          X=[M,nei_prev_bug_basic(:,g), ego_prev_bug,ones(filenum*(release_num-1),1)];
          index= index_test & (ego_curr_bug>=0) & (nei_prev_bug_basic(:,g)<6);
      end
      
      varnames={'number of changes','number of authors','Cyclomatic','Essent_Complex','LOC','prev_nei_bug','prev_ego_bug','constant'};
      [betahat, alphahat, results] = gee(id(index), y(index), t(index), X(index,:), 'b','equi',varnames);
      GEE=[results.model{8,2},results.model{8,6},results.model{8,7},results.model{8,5}];
     
    
    end

     
   end
