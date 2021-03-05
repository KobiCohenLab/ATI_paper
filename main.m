

% Simulations for active sensing of epidemics with delyed detections
% This is a stochastic model in contrast to SIR models
% Testing detection and quarantining strategies
clear all
%load_initialization_flag=0;

nj=107
vsampling=1:5;
ntest=length(vsampling);
save_flag=1;
movie_flag=0;
 p_long_contact=0.5;   
forget_factor_coeff=0.7; % forgetting factor on belief between 0 to 1, 0 - no memory between time steps, 1- full memory
forget_factor_test=0.25; %for test_rate=0.01, p_detect=1, factor 0.2 was the best, for p_detect 0.8, 0.4 was the best
    
    % T=1000;
T=500;
trials=5;
r_infect=2.3; % number of people in average that infected person can infect
N0=50; %initial number of infected people
p_asymptomatic=0.3;%.3; %15; %probability of an infected person to be asymptomatic
    p_severe=0.05; %probability of an infected person to be in severe state
    p_death=0.3;%probability of severe infected people to die
    
    graph_type=1;% Random Geometric graphs
    load_graph_flag=0;
    N=1e5
    if load_graph_flag
        load graph_random_geometric;
    else
        %    Generate_graph(graph_type)
R=40000;



location=rand(N,2)*R;
  %      save graph_random_geometric
    end
    
    % Generate distance map
    kk=1;
    p_hub=0.05; %probability of user to be able to have far neighbors (used to geerate Drand)
    
    test_rate=0.03; %percents of the population that we test
    vN=[200 100 50 25];
    %vN=[25 25 25 25 ];
    vsig=[40 20  10  5];
    vp_degree=[0 0.01 0.03 0.15 1];
    
    lvN=length(vN);
    lvp_degree=length(vp_degree)
    rand_degree=rand(N,1);
    v_degree=zeros(N,1);
    for mm=1:lvp_degree-1
        ind=find((rand_degree<vp_degree(mm+1)) & (rand_degree>vp_degree(mm)));
        v_degree(ind)=vN(mm);
    end
    
    
    v_hub=(rand(N,1)<p_hub); %N x 1 vector of far neighbors (1 indicates hub, zero indicates no) (used to geerate Drand)
    
    %make everything symmetry: p_H0_H0=p_H1_H1=p_detect
    p_detect=0.8;
    %p_detect=1;
    
  
    
    Drand=R*ones(N,1);
    
    %Drand=2e3*(1-v_neighbors).*ones(N,1)+R*v_neighbors.*ones(N,1); %N x 1 distance vector . Drand(ii) is the radius for agent ii
    %(Drand(ii)=R for user ii that has far neighbors, Drand(jj)=2e3 for user jj that has only near neighbores)
    %Drand=R*ones(N,1);;
    
 
    
    
    DeltaT_min_infection_symptoms=2;%minimal time from infection until symptoms appear (for symptomatic people)
    DeltaT_max_infection_symptoms=11;%maximal time from infection until symptoms appear (for symptomatic people)
    DeltaT_symptoms_quarantine=3;
    
    DeltaT_min_infection_infectious=1;%minimal time from infection until ainfected people start to infect
    
    DeltaT_quar=14;
    Tdetection_delay=1;
    %v_p_quar_success=[0:0.2:1]; %probability that a person in qurantine will not infect
    v_p_quar_success=[0.7 0.8 0.9 1]
    lvp=length(v_p_quar_success);
    
    %v_p_quar_success=[1]
    total_infected=zeros(T+1,trials,lvp,ntest);
    total_heal=zeros(T+1,trials,lvp,ntest);
    total_quar=zeros(T+1,trials,lvp,ntest);
    total_test=zeros(T+1,trials,lvp,ntest);
    v_heal=zeros(N,1);
    v_positie=zeros(N,2);
    max_neighbors=600;
    for i_p_quar_success=1:lvp %plotting as a function of p_success
        p_quar_success=v_p_quar_success(i_p_quar_success);

        for tr=1:trials %averaging over trails
            tr
            if tr==20*floor(tr/20)
                tr
            end
            
            for ii=1:N
                nu=max(1,min(max_neighbors,round(v_degree(ii)+randn(1)*(v_degree(ii)/7))));
                %        ik=find(iu>=ii);
                %        iu(ik)=iu(ik)+1;
                
                A(ii).ngb=randperm(N,nu);
                A(ii).l_ngb=nu;
            end
            for ii=1:N
                for jj=1:A(ii).l_ngb
                    kk=A(ii).ngb(jj);
                    if ~sum(ii==A(kk).ngb);
                        lA=A(kk).l_ngb;
                        A(kk).l_ngb=lA+1;
                        A(kk).ngb(lA+1)=ii;
                    end
                    
                    
                end
                
            end
            initial_ind=randperm(N); 
         itest=0;
        for sampling_flag= vsampling; %[1 2 3 4];% [1 2 3] 
        itest=itest+1;
          detect2_flag=0;
    if sampling_flag==1
        detect2_flag=1;
    end
      if sampling_flag==4
        p_H0_H0=1;
        p_H1_H1=0;
    else
        p_H1_H1=p_detect; %declared infected given infected
        p_H0_H0=1;%-0.01; %declared healthy given healthy
    end
      p_H1_H0=1-p_H0_H0; %declared infected given healthy (false positive)
      p_H0_H1=1-p_H1_H1; %declared healthy given infected
  
        % clearvars -except T trials v_p_quar_success v_avg_l_infected mat_avg_l_severe mat_avg_l_test v_avg_l_quar p_quar_success
       
            for ii=1:N
                
                Neffective_neighbors(ii,1)=A(ii).l_ngb;
                p_infection(ii)=r_infect/Neffective_neighbors(ii); %probability of infection to yield r_infect in average.
            end % ii                                                %should be dependent on Mean_Neffective_neighbors in the future
            %Generate agents
            vT_heal=zeros(N,1);
            for ii=1:N;
                %     if ii==1000*floor(ii/1000) %stopping condition for debugging
                %         ii
                %     end
                
                % graph_structure
                
                lngb=A(ii).l_ngb;
                agent(ii).l_ngb=lngb;
                agent(ii).ngb=A(ii).ngb;
                agent(ii).pngb=rand(lngb,1)<p_long_contact;
                nn=location(ii,1);%coordinate x of agent i
                mm=location(ii,2);%coordinate y of agent i
                agent(ii).x=nn;%coordinate x of agent i
                agent(ii).y=mm;%coordinate y of agent i
                
                %kk=find(sqrt((location(:,1)-nn).^2+(location(:,2)-mm).^2)<=Drand(ii));%#neighbors x 1 vector, find all agent ii's neighbors in distance Drand(ii)
                %s=randperm(length(kk));
                %                 if length(kk)<Neffective_neighbors(ii) %agent(ii).ngb chooses only 100 neighbors if the random number of neighbors exceeds Neffective_neighbors
                %                     agent(ii).ngb=kk;
                %                 else
                %                     agent(ii).ngb=kk(s(1:Neffective_neighbors(ii)),1);
                %                 end
                %                 % state variables
                
                
                
                %agent(ii).l_ngb=length(agent(ii).ngb); %number of neighbors
                agent(ii).infected=0; %0=infected, 1=healthy
                agent(ii).detected=0; %0=not-detected, 1=detected
                agent(ii).quarantine=0; %0=not-quarantine, 1=quarantine
                agent(ii).heal=0; %0=not-quarantine, 1=quarantine
                
                agent(ii).infectious=0; %0=not-infecious, 1=infecious
                agent(ii).asymptomatic=rand(1)<p_asymptomatic; %0=not-asymptmotic, 1=asymptomatic
                agent(ii).severe=0; %0=not severe, 1=severe
                %%uniform dist of symptoms time:
                %agent(ii).dt_symptoms=round(DeltaT_min_infection_symptoms+rand(1)*(DeltaT_max_infection_symptoms-DeltaT_min_infection_symptoms));
                %Rayliegh dist of symptoms time:
                agent(ii).dt_symptoms=round(raylrnd(5.1/sqrt(2*log(2))));
                if agent(ii).dt_symptoms>20
                    agent(ii).dt_symptoms=20;
                end
                %random time where symptoms appear (between DeltaT_min_infection_symptoms to DeltaT_max_infection_symptoms)
                agent(ii).dt_infectious=14; %total time the agent is infectious
                %agent(ii).dt_symptoms=3;
                % time variables
                agent(ii).t_infect=0; % absolute infection time
                agent(ii).t_detect=0; % absolute detection time
                agent(ii).t_heal=0;   % absolute detection time
                agent(ii).t_death=0;  % absolute detection time
                agent(ii).ngb_quarantine=0;
                agent(ii).qtested=0;
                agent(ii).positive=0;
                 vT_heal(ii)=21+round(10*rand(1));
                 agent(ii).DeltaT_heal=vT_heal(ii);
            end
            
            % Initialize
            v_infected=zeros(N,2); % N x 2 vector of infected people, First coordinae is identity (1=sick, 0=healthy).
            %                                  Second coordinate is infection time
            v_detected=zeros(N,2); % N x 2 vector of detected people, First coordinae is identity (1=detected, 0=not-detected).
            %                                  Second coordinate is detection time
            v_quar=zeros(N,2); % N x 2 vector of people in quarantine, First coordinae is identity (1=quar, 0=not).
            %                                  Second coordinate is time they were sent to qurantine
            v_heal=zeros(N,2); % N x 2 vector of healed people, First coordinae is identity (1=healed, 0=not).
            %                                  Second coordinate is time they were sent to qurantine
            real_quar=zeros(N,1);
            
            list_of_infected_neighbors=zeros(N,1);
            kk=0;
                       ind=initial_ind;
            for jj=1:N0; %initilizing first N0 sick agents
                ii=ind(jj);
                agent(ii).infected=1;
                agent(ii).t_infect=0;; %infection time for agent ii, first N0 agents are infected at time t=0
                
                v_infected(ii,1)=1; %assigning infection identity to v_infected
                v_infected(ii,2)=0; %assigning infection time to v_infected
                kk=kk+1;
                list_of_infected_nodes(kk)=ii; %first it contains only first N0 agents
            end
            l_infected=sum(v_infected(:,1)); %number of infected agents, need only first coordinate (second coordinate is sum of infection time)
            
            %save graph_random_geometric
            N_sympt=0;
            v_test_belief=zeros(N,1);
            v_positive=zeros(N,2);
            lvd=T;
            F(lvd) = struct('cdata',[],'colormap',[]);
            for t=1:T
                tic
                if t==60
                    1;
                end
                %v_test_belief=zeros(N,1); % N x 1 belief vector, for each agent, counts the numer of its 2nd cycles
                v_test_belief=forget_factor_coeff*v_test_belief; %adding forgetting factor to belief
                % Move healed agents to healed agents list
                list_of_infected_nodes=find(v_infected(:,1));
                l_infected=sum(v_infected(:,1));
                for nn=1:l_infected
                    nagent=list_of_infected_nodes(nn);
                    if t==agent(nagent).t_infect+vT_heal(nagent);
                        agent(nagent).heal=1;
                        agent(nagent).infected=0;
                        agent(nagent).quarantine=0;
                        real_quar(nagent)=0;
                        agent(nagent).t_heal=t;
                        v_heal(nagent,1)=1;
                        v_infected(nagent,1)=0;
                        v_quar(nagent,1)=0;
                        v_quar(nagent,2)=0;
                        real_quar(nagent)=0;
                    end
                    % release from quarantine
                end
                
                for nagent=1:N
                    if t==v_quar(nagent,2)+DeltaT_quar && ~agent(nagent).infected;
                        agent(nagent).qtested=0;
                        v_quar(nagent,1)=0;
                        v_quar(nagent,2)=0;
                        agent(nagent).quarantine=0;
                        real_quar(nagent)=0;
                        v_test_belief(nagent)=0;
                        
                    end
                end
              ;
                %l_quar=sum(v_quar(:,1));
                l_quar=sum(real_quar);
%                 for ii=1:N
%                     nq=agent(ii).quarantine;
%                        l_quar=l_quar+nq;
%                  end
                l_heal=sum(v_heal(:,1));
                total_infected(t,tr,i_p_quar_success,itest)=l_infected;
                total_quar(t,tr,i_p_quar_success,itest)=l_quar;
                total_heal(t,tr,i_p_quar_success,itest)=l_heal;
                
                for nn=1:l_infected
                    nagent=list_of_infected_nodes(nn); %nagent=identity number of agent
                    if ~agent(nagent).quarantine % if not in quarantine it can infect other agents
                        % Set infected nodes and infection times by a newly infected agent.
                        %  if ~isempty(list_of_infected(nn))
                        if t-1==agent(nagent).t_infect % then we randomize the infection times of all infected neighbors
                            vp_infection=p_infection(nagent).*(0.25+1.75*agent(nagent).pngb);
                            v_ind=(rand(agent(nagent).l_ngb,1)<vp_infection*(min(1,agent(nagent).dt_symptoms/4)));
                           
                            
                            %#(agent's nagent neighbors) x 1 vector, 1=infected, 0=not infected
                            %*** issue: the effective number of neighbors is less than 100, we might want to
                            %slightly increase p_infection
                            agent(nagent).vector_of_infected_neighbors=v_ind;
                            ind1=find(v_ind==1);
                            %number of infected neighbors:
                            agent(nagent).number_of_infected_neighbors=sum(v_ind);
                            %list of infected neighbors of agent nagent:
                            agent(nagent).list_of_infected_neighbors=agent(nagent).ngb(ind1);
                            l1=agent(nagent).number_of_infected_neighbors;
                            %updating the random times where neigbors are infected
                            %                            agent(nagent).list_of_infection_times=t+round(DeltaT_min_infection_infectious+rand(l1,1)*(agent(nagent).dt_symptoms+DeltaT_symptoms_quarantine-DeltaT_min_infection_infectious));
                            agent(nagent).list_of_infection_times=t+round(DeltaT_min_infection_infectious+rand(l1,1)*((agent(nagent).dt_symptoms+DeltaT_max_infection_symptoms)/2-DeltaT_min_infection_infectious));
                            %                            agent(nagent).list_of_infection_times=t+round(DeltaT_min_infection_infectious+rand(l1,1)*(agent(nagent).dt_symptoms+1-DeltaT_min_infection_infectious));
                            %                            agent(nagent).list_of_infection_times=t+round(DeltaT_min_infection_infectious+rand(l1,1)*(DeltaT_max_infection_symptoms-DeltaT_min_infection_infectious));
                            
                        end %end if
                        
                        symptoms_time=agent(nagent).t_infect+agent(nagent).dt_symptoms;
                        
                        if (t<=symptoms_time+DeltaT_symptoms_quarantine || agent(nagent).asymptomatic)%  && ~agent(nagent).positive % update infections only if not yet detected
                            %updating the infecected neighbors of agent nagent:
                            for kk=1:agent(nagent).number_of_infected_neighbors
                                current_infected=agent(nagent).list_of_infected_neighbors(kk);
                                agent_heal_flag=agent(current_infected).heal;
                                agent_quarantine_flag=agent(current_infected).quarantine;
                                agent_infected_flag=agent(current_infected).infected;
                                if agent(nagent).list_of_infection_times(kk)==t && ~agent_heal_flag && ~agent_infected_flag && ~agent_quarantine_flag
                                    agent(current_infected).t_infect=t;
                                    agent(current_infected).infected=1;
                                    v_infected(current_infected,1)=1;
                                    v_infected(current_infected,2)=t;
                                    
                                    
                                end
                            end
                        end
                    end
                end
    
                list_of_update=find(v_infected(:,1) | v_positive(:,1));
                l_update=length(list_of_update);
              
                for nn=1:l_update
                    nagent=list_of_update(nn);
                    infected_flag=agent(nagent).infected;
                    symptoms_time=agent(nagent).t_infect+agent(nagent).dt_symptoms;
                    symptoms_flag=(t==symptoms_time+DeltaT_symptoms_quarantine);
                    detection_time_flag=( t==v_positive(nagent,2)+Tdetection_delay);
                    positive_flag=(v_positive(nagent,1) && detection_time_flag);
                    if  positive_flag || (infected_flag && symptoms_flag)
                            if ~v_quar(nagent,1);
                               qrand=(rand(1)<p_quar_success);
                                agent(nagent).quarantine=qrand;
                                real_quar(nagent)=qrand;
                               
                                v_quar(nagent,1)=1;
                                v_quar(nagent,2)=t;
                            end        
                         
                        if ~agent(nagent).ngb_quarantine 
                            agent(nagent).ngb_quarantine=1;
                            lneighbors=agent(nagent).l_ngb;
                            for kk=1:lneighbors
                                i_quar=agent(nagent).ngb(kk);%index of neighbore in quarantine
                                contact_type=agent(nagent).pngb(kk);
                                p_specific=0.5+contact_type;
                                if ~v_quar(i_quar,1) && ~agent(i_quar).heal
                                    if sampling_flag<4
                                    v_test_belief(i_quar,1)=min(1,v_test_belief(i_quar,1)+r_infect/lneighbors*p_specific);
                                    end
                                    if sampling_flag==5
                                      v_test_belief(i_quar,1)=1; %min(1,v_test_belief(i_quar,1)+r_infect/agent(nagent).l_ngb);
                                    end
                                      qrand=(rand(1)<p_quar_success);
                                    agent(i_quar).quarantine=qrand; %neighbor in quarantine with probability p_quar_success
                                    real_quar(i_quar,1)=qrand;
                                    v_quar(i_quar,1)=1; %vector of all agent in quaranteine and times the y sent to quarantive
                                    v_quar(i_quar,2)=t;
                                     
                                end
                                
                                % update belief on all its neighbors
                                if detect2_flag
                                    lneighbors1=agent(i_quar).l_ngb;
                                    tot_neighbors=lneighbors*lneighbors1;
                                    for mm=1:lneighbors1
                                        i_detect= agent(i_quar).ngb(mm);
                                         contact_type1=agent(i_quar).pngb(mm);
                                         p_specific1=0.5+contact_type1;
                                        if ~v_quar(i_detect,1) && ~agent(i_detect).heal
                                            %            v_test_belief(i_detect,1)=v_test_belief(i_detect,1)+1; %we update the 2nd circle count i_detect
                                            v_test_belief(i_detect,1)=min(1,v_test_belief(i_detect,1)+r_infect^2/( tot_neighbors)*p_specific*p_specific1); %we update the 2nd circle count i_detect
                                           

                                        end
                                    end
                                end
                                
                            end
                            

                        end
                    
                    end 
                    
                    
                end
                
                
                %Testing all agenets with high belief----
if sampling_flag==1 || sampling_flag==2 || sampling_flag==3 || sampling_flag==5
                    No_test=round(test_rate*N); %number of tests at time t
                end
                if sampling_flag==1 || sampling_flag==3 || sampling_flag==5
                    [v_test_belief_sort  v_test_belief_sort_i]=sort(v_test_belief,'descend');
                end
                if sampling_flag==2 %sampling_flag=3 means no detection, just need to declare v_test_belief_sort_i
                    v_test_belief_sort_i=randperm(length(v_test_belief));
                end
                
                test_count=1;
                break_count=1;
              
                    
                    
                
                
                if sampling_flag==1  || sampling_flag==3 ||   sampling_flag==2 ||   sampling_flag==5
                    while test_count<=No_test && break_count<=N
                        n_agent=v_test_belief_sort_i(break_count);
                        break_count=break_count+1;
                        flag_1=~agent(n_agent).heal;
                        flag_2=~v_quar(n_agent,1);
                        flag_3= ~agent(n_agent).qtested;
                         if  flag_1 &&  (flag_2 || flag_3);
                        %    if  flag_1 &&  (flag_2 )
                              total_test(t,tr,i_p_quar_success)=total_test(t,tr,i_p_quar_success)+1;
                            test_count=test_count+1;
                            if v_quar(n_agent,1);
                                agent(n_agent).qtested=1;
                            end
                            if agent(n_agent).infected
                                if rand(1)<p_H1_H1 %declare H1
                                   if ~v_quar(n_agent,1)
                                       qrand=rand(1)<p_quar_success;
                                       agent(n_agent).quarantine=qrand;
                                       real_quar(n_agent)=qrand;
                                   end
                                    agent(n_agent).positive=1;
                                    
                                    v_positive(n_agent,1)=1;
                                    v_positive(n_agent,2)=t;
                                    v_quar(n_agent,1)=1; %vector of all agent in quaranteine and times the y sent to quarantive
                                    v_quar(n_agent,2)=t;
                                    
                                else %declare H0
                                    v_test_belief(n_agent)=forget_factor_test*v_test_belief(n_agent); %adding forgetting factor to belief
                                    
                                end
                            else
                                if rand(1)>p_H0_H0 %declare H1
                                    if ~v_quar(n_agent,1)
                                        qrand=rand(1)<p_quar_success;
                                        real_quar=qrand; 
                                        agent(n_agent).quarantine=qrand;
                                         
                                    end
                                         agent(n_agent).positive=1;
                                    v_positive(n_agent,1)=1;
                                    v_positive(n_agent,2)=t;
                                    v_quar(n_agent,1)=1; %vector of all agent in quaranteine and times the y sent to quarantive
                                    v_quar(n_agent,2)=t;
                                   
                                else %declare H0
                                    v_test_belief(n_agent)=forget_factor_test*v_test_belief(n_agent); %adding forgetting factor to belief
                                end
                            end
                            
                            
                        end
                    end
                end
                
                
                
                
                
                
                
                %-----------------------------------------
                l_infected=sum(v_infected(:,1)); %==total_infected(t+1)
                l_quar=sum(v_quar(:,1));
                
                list_of_infected_nodes=find(v_infected(:,1));
                if movie_flag
                    figure(1)
                    hold off
                    plot(location(:,1).*(1-(v_quar(:,1)|v_infected)),location(:,2).*(1-(v_quar(:,1)|v_infected)),'go');
                    title(['t=' int2str(t)]);
                    hold on
                    plot(location(:,1).*v_quar(:,1),location(:,2).*v_quar(:,1),'b.');
                    plot(location(:,1).*v_infected(:,1),location(:,2).*v_infected(:,1),'k*');
                    plot(0,0,'go');
                    
                    %pause(0.1)
                    drawnow
                    
                    F(t) = getframe;
                    hold off
                    
                end
                
                
            end %end for t=1:T
            if movie_flag
                writerObj = VideoWriter(['Fsamp' int2str(sampling_flag) '.mp4'],'MPEG-4')
                writerObj.Quality=100;
                open(writerObj)
                writeVideo(writerObj,F);
                close(writerObj)
                %                movie2avi(F,'Spread.avi','Compression','None')
            end
            
            
            %         total_severe=total_infected*p_severe; %T x 1 vextor of sever infected agents with time
            %         mat_l_severe(:,tr)=total_infected*p_severe; %T x 1 vextor of sever infected agents with time
            %         v_l_infected(tr)=l_infected %==total_infected(t+1)
            %         v_l_quar(tr)=l_quar
            sampling_flag
             p_quar_success
            end %for sampling_flag
        if tr==10*floor(tr/10)
              if save_flag
                  %eval(['save data_all_samples_real_belief_j' int2str(nj)  'T' int2str(test_rate*1000) 'Q' int2str(100*p_quar_success) ]); 
                      save_file_name_str=['save data_all_samples_belief1_j' int2str(nj)  'T' int2str(round(test_rate*1000)) 'Q' int2str(100*p_quar_success)  'N' int2str(round(N/1000)) 'Ni'  int2str(N0)]
                     eval(save_file_name_str);
              end
    end
            
        end %for trials
        %     v_avg_l_severe=mean(mat_l_severe,2); %avergaing severe over trials to get vector of #sever vs time
        %     mat_avg_l_severe(find(v_p_quar_success==p_quar_success),:)=v_avg_l_severe';
        %     %how many infected on average:
        %     avg_l_infected=mean(v_l_infected);
        %     v_avg_l_infected(find(v_p_quar_success==p_quar_success))=avg_l_infected
        %
        %     %how many tests have been done on average
        %     %avg_l_test=mean(v_l_test);
        %     %v_avg_l_test(find(v_p_quar_success==p_quar_success))=avg_l_test
        %
        %     %how many people in quartine on average
        %     avg_l_quar=mean(v_l_quar);
        %     v_avg_l_quar(find(v_p_quar_success==p_quar_success))=avg_l_quar
        
    end %for v_p_quar_success
    % mat_avg_l_severe;
    % v_avg_l_infected %how many infected on average
    % v_avg_l_quar %how many people in quartine on average
    
    %save vs_p_success_var mat_avg_l_severe v_avg_l_infected v_p_quar_success p_hub T%saving vs_p_success variables
    %plotting;
    %grid on
    
     
       M_infected=squeeze(total_infected(:,:,1,:));
       M_quar=squeeze(total_quar(:,:,1,:));          
       M_heal=squeeze(total_heal(:,:,1,:));
     sort_M_infected=sort(M_infected,2,'descend');
   sort_M_quar=sort(M_quar,2,'descend');
   sort_M_heal=sort(M_heal,2,'descend');
   
     figure(1)
     vt=1:(T+1);
          plot(vt,squeeze(mean(sort_M_infected(:,:,1),2)),'b-'); 
          hold on
         
          plot(vt,squeeze(mean(sort_M_infected(:,:,2),2)),'r--');
          plot(vt,squeeze(mean(sort_M_infected(:,:,3),2)),'g.-');
          plot(vt,squeeze(mean(sort_M_infected(:,:,4),2)),'k+-');
             plot(vt,squeeze(mean(sort_M_infected(:,:,5),2)),'m.-');
          legend('Active 2', 'random','Active 1', 'no detection','DTI')
hold off
     figure(2)
             plot(vt,squeeze(mean(sort_M_quar(:,:,1),2)),'b-'); 
          hold on
        
          plot(vt,squeeze(mean(sort_M_quar(:,:,2),2)),'r--');
          plot(vt,squeeze(mean(sort_M_quar(:,:,3),2)),'g.-');
          plot(vt,squeeze(mean(sort_M_quar(:,:,4),2)),'k+-');
          plot(vt,squeeze(mean(sort_M_quar(:,:,5),2)),'m.-');
          legend('Active 2', 'random','Active 1', 'no detection','DTI')
hold off
          figure(3)
              plot(vt,squeeze(mean(sort_M_heal(:,:,1),2)),'b-'); 
          hold on
          
          plot(vt,squeeze(mean(sort_M_heal(:,:,2),2)),'r--');
          plot(vt,squeeze(mean(sort_M_heal(:,:,3),2)),'g.-');
          plot(vt,squeeze(mean(sort_M_heal(:,:,4),2)),'k+-');
          plot(vt,squeeze(mean(sort_M_heal(:,:,5),2)),'m.-');
          legend('Active 2', 'random','Active 1', 'no detection','DTI')
hold off
    if save_flag
        save_file_name_str=['save data_all_samples_belief1_j' int2str(nj)  'T' int2str(round(test_rate*1000)) 'Q' int2str(100*p_quar_success)  'N' int2str(round(N/1000)) 'Ni'  int2str(N0)]
         eval(save_file_name_str);
      
    end
    


avg_total_infected=mean(total_infected,2);
avg_total_quar=mean(total_quar,2);
avg_total_heal=mean(total_heal,2);

% figure(1)
% plot(1:T+1, avg_total_infected(:,1,1)/N, 'b', 1:T+1, avg_total_infected(:,1,2)/N, 'k', 1:T+1, avg_total_infected(:,1,3)/N, 'r' )
% grid on;
% title('Total Infected Rate')
% legend('quarantine success rate=0.6', 'quarantine success rate=0.8', 'quarantine success rate=1')
%
% figure(2)
% plot(1:T+1, avg_total_heal(:,1,1)/N, 'b', 1:T+1, avg_total_heal(:,1,2)/N, 'k', 1:T+1, avg_total_heal(:,1,3)/N, 'r' )
% grid on;
% title('Total Heal Rate')
% legend('quarantine success rate=0.6', 'quarantine success rate=0.8', 'quarantine success rate=1')
%
% figure(3)
% plot(1:T+1, avg_total_quar(:,1,1)/N, 'b', 1:T+1, avg_total_quar(:,1,2)/N, 'k', 1:T+1, avg_total_quar(:,1,3)/N, 'r' )
% grid on;
% title('Total Qaurantine Rate')
% legend('quarantine success rate=0.6', 'quarantine success rate=0.8', 'quarantine success rate=1')
