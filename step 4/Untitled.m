%% /************************************************************************************************/
%  /*** Topic: Sync calculation                                                                  ***/
%  /***                                                                         Ali-Seif         ***/
%  /*** Version Release 9.8.0.1323502                                                            ***/
%  /*** Date: 11/29/2020                                                                         ***/
%  /*** Code implemented in MATLAB R2020a 64bit                                                  ***/
%  /*** MSI: PX60 6QD/ DDR4                                                                      ***/
%  /*** Run under a Intel® Core™ i7-6700HQ CPU @ 2.60GHz × 64 based processor with 16 GB RAM     ***/
%  /************************************************************************************************/
%% ##############################################################
%  ####                                                      ####
%  ####                   Cleaning                           ####
%  ####                                                      ####
%  ############################################################## 
clc;                                                            %Clear Command Window
clear;                                                          %Remove items from workspace,freeing up system memory
clf;                                                            %Clear current figure window
%% ##############################################################
%  ####                                                      ####
%  ####            Raed file and convert to array            ####
%  ####                                                      ####
%  ##############################################################
N=100;                                                            %number of neurons
phi=1;
I=2;
File=dlmread('C:\Users\Ali\Desktop\temp.txt');                                       %raed file
t=File(:,1);                                                    %convert file to time
v = zeros(length(File),N);                                      %craet fill matrix with length File and N
for  i=1:1:N
    v(:,i)=File(:,i+1);                                         %convert file to voltage
end
s_total=File(:,N+2);                                            %convert file to total wight synaptic
%% ##############################################################
%  ####                                                      ####
%  ####          plot voltage .VS time for any neurons       ####
%  ####                                                      ####
%  ##############################################################
tiledlayout(4,1)                                                %creat 4 figure in one page
nexttile                                                        %go to the first figure
for  i=1:1:N
    plot(t,v(:,i))                                              %plote who to each neurons spike during the time
    hold on                                                     %plot the next figure on this figure
end
title (['Number of neurons(N)=',num2str(N),'Possibility of connection(P)=',num2str(0.3),'External Current(I)=',num2str(I)]);
ylabel( 'Voltage(mv)' );
set(gca,'xlim',[min(t) max(t)],'ylim',[-85 60])                 %Figure framework
fid=fopen(['C:\Users\Ali\Desktop\spike_',num2str(phi),'_',num2str(I),'.txt'],'w');        
                                                                %create File .txt 
nexttile                                                        %go to the second figure
%% ##############################################################
%  ####                                                      ####
%  ####    find exact time fier point and plot rater plot    ####
%  ####                                                      ####
%  ##############################################################
y2 = sin(t)-sin(t)-52.0;                                        %Select the threshold line
c = ones(1,0);                                                  %create fill array
for number_neuron=1:1:N
    [xout,yout] = intersections(t,v(:,number_neuron),t,y2,1);   %when voltage Curve for each neurons cross threshold Report x and y point
    if rem(length(xout),2) == 0           
            n=0.5*(length(xout));                               %if number of point report was even 
        else           
            n=0.5*(length(xout)+1);                             %if number of point report was odd
    end
    cxout = zeros(1,n);                                         %craet fill array with length n
    for i=1:1:n
             cxout(1,i)=xout(2*i-1);                            %Save points one by two
    end
     a=repmat (cxout(1,:),2 );                                  %two point for each spike point
     b=[number_neuron-0.5,number_neuron+0.5];                   %length Vertical rater plot for each neuron
     plot (a, b, "k-",'linewidth',2);                           %plot rater plot for each neuron
     c = [c cxout(1,:)];                                        %add this array to before array
    hold on                                                     %plot the next figure on this figure
end
ylabel( 'Node index' );
set(gca,'xlim',[min(t) max(t)],'ylim',[0 number_neuron+1])      %Figure framework
nexttile                                                        %go to the next figure
%% ##############################################################
%  ####                                                      ####
%  ####plot rater plot for all neurons and sum wight synaptic####
%  ####                                                      ####
%  ##############################################################
a=repmat (c,2 );                                                %two point for each spike point of all neurons
b=[0.5,1.5];                                                    %length Vertical rater plot for all spike line
plot (a, b, "k-",'linewidth',0.5);                              %plot rater plot for all neuron
ylabel( 'Spike train' );
set(gca,'xlim',[min(t) max(t)],'ylim',[0 2])                    %Figure framework
nexttile                                                        %go to the next figure
plot(t,s_total)                                                 %plot sum wight synaptic during the time
xlabel( 'Time(ms)' );
ylabel( 'Sygma(S_i)' );
%% ##############################################################
%  ####                                                      ####
%  ####            Rad file and convert to array             ####
%  ####                                                      ####
%  ##############################################################
c=transpose(c);                                                 %transpose matrix
c=sortrows(c);                                                  %Arrange the spikes in order low to high
c=transpose(c);                                                 %transpose matrix
sort = [length(c) N phi I];                                          %creat array with 2 paramiter
sort = [sort c];                                                %add point spike array to before array
sort=transpose(sort);                                           %transpose matrix
fprintf(fid, '%f\n', sort');                                    %print in to the file
fclose(fid);                                                    %close file .txt
N
