#This script compares pdb files and write pairwise RMSD in a table: total_structs
#Here it's assumed that the name of the pdb files are struct# with #=00001,00002 and so on.
#It's also assumed that the pdb files are written in a certain way like in
#the ensemble with name struct#. Functions readpdb and readpdb1 might need to have
#inputs adjusted if the pdb file is a little bit different. Like in
#ensemble mnodes. 

#Note: File is not yet complete


pdbfilenameA = 'G:\My Drive\University of Maryland\PhD Project\Ensembles\STRUCT_superimposed_A\struct00001.pdb'
pdbfilenameB = 'G:\My Drive\University of Maryland\PhD Project\Ensembles\STRUCT_superimposed_A\struct00001.pdb'
#pdbfilenameA = 'G:\My Drive\University of Maryland\Summer Project 2020\Ensembles\MNODES_superimposed_B\mnodesS00001.pdb'
#pdbfilenameB = 'G:\My Drive\University of Maryland\Summer Project 2020\Ensembles\MNODES_superimposed_B\mnodesS00001.pdb'

struct = strfind(pdbfilenameB,'struct00001') + 10#FINDS THE POSITION OF 1 IN THE FILE NAME IN ORDER TO BE ABLE TO CHANGE FOR 2,3,4 ...
#struct = strfind(pdbfilenameB,'structS00001') + 11 
#struct = strfind(pdbfilenameB,'mnodesS00001') + 11

numstruct = length(keep68) #TOTAL NUMER OF STRCTURES IT'S GOING TO BE USED

off = 1 #for STRUCT, 0 for MNODES
#read in first domain
[coorA1,atnamA1,at_resA1]=readpdb(pdbfilenameA,[],[],1,off,'A')
#read in second domain
[coorA2,atnamA2,at_resA2]=readpdb(pdbfilenameA,[],[],1,off,'B')

[coorB1,atnamB1,at_resB1]=readpdb(pdbfilenameB,[],[],1,off,'A')
[coorB2,atnamB2,at_resB2]=readpdb(pdbfilenameB,[],[],1,off,'B')

total_structs = zeros(length(keep68),length(keep68))

for m in range (numstruct):
    i = keep68(m)
    pdbfilenameA = 'G:\My Drive\University of Maryland\PhD Project\Ensembles\STRUCT_superimposed_A\struct00001.pdb'
    if (i < 10):
        pdbfilenameA(struct) = num2str(i)
    elif (i < 100):
        pdbfilenameA(struct-1:struct) = num2str(i)
        else
            if (i < 1000):
                pdbfilenameA(struct-2:struct) = num2str(i)
            else
                if (i < 10000):
                    pdbfilenameA(struct-3:struct) = num2str(i)
                else
                    if (i < 100000):
                        pdbfilenameA(struct-4:struct) = num2str(i)
                        
    coorA1=readpdb1(pdbfilenameA,[],[],1,at_resA1,at_resA1) #FOR STRUCTS
    coorA2=readpdb1(pdbfilenameA,[],[],1,at_resA2,at_resA2) #FOR STRUCTS
    #[coorA1,atnamA1,at_resA1]=readpdb(pdbfilenameA,[],[],1,off,'A') #FOR MNODES
    #[coorA2,atnamA2,at_resA2]=readpdb(pdbfilenameA,[],[],1,off,'B') #FOR MNODES
    
    a = 1
    for s = 1:numstruct
        j = keep68(s)
        #total_structs(m,1) = i
        
        if j ~= i

            if j < 10
                pdbfilenameB(struct-4:struct) = "0000" + num2str(j)
            else
                if j < 100
                    pdbfilenameB(struct-4:struct) = "000" + num2str(j)
                else
                    if j < 1000
                        pdbfilenameB(struct-4:struct) = "00" + num2str(j)
                    else
                        if j < 10000
                            pdbfilenameB(struct-4:struct) = "0" + num2str(j)
                        else
                            if j < 100000
                                pdbfilenameB(struct-4:struct) = num2str(j)
                            end
                        end    
                    end           
                end
            end    

            #The other picture
            #coorB1=readpdb1(pdbfilenameB,[],[],1,at_resB1,at_resB1) #FOR STRUCTS
            #coorB2=readpdb1(pdbfilenameB,[],[],1,at_resB2,at_resB2)
            [coorB1,atnamB1,at_resB1]=readpdb(pdbfilenameB,[],[],1,off,'A')#FOR MNODES
            [coorB2,atnamB2,at_resB2]=readpdb(pdbfilenameB,[],[],1,off,'B')

            #list of residues to superimpose
            reslst=[]

            [coorBrot,Rotmat,rmsd,coor2s]=superimpose([coorA1coorA2],[atnamA1atnamA2],[at_resA1at_resA2],[coorB1coorB2],[atnamB1atnamB2],[at_resB1at_resB2],reslst,reslst,[],[coorB1coorB2])

            total_structs(m,a) = rmsd
            a = a+1
        end
    end
    #save total_structs_rmsd_STRUCT_Clust_2A_maxweight.mat total_structs '-v7.3'
end

#This loop is just to put the vector total_structs that has the pairwise
#RMSD as a matrix with zeros in the diagonal. 
for i = 1:length(total_structs(:,1))
    a = nonzeros(total_structs(i,:))
    total_structs(i,i) = 0
    total_structs(i,i+1:end) = a(i:end)
end

[rows columns] = size(total_structs) 
v = [] #v is going to be the vector that contains all the pairwise RMSD
for i = 1:rows-1 
    v = [v total_structs(i,i+1:columns)] 
end

tree=linkage(v,'average') #ths function linkage does hierarchal clustering

dendrogram(tree,0) #this plots the tree with the clusterings

T = cluster(tree,'cutoff',4,'Criterion','distance') #the RMSD clustering is chose to be 4A in here
T(:,2) = 1:length(T) 
T(:,3) = keep68
T = sortrows(T)

xsol = result_of_maxent.x(:,index68)

#Here I choose the representative of the cluster. The way I do, the
#representative is the structure with the highest initial weight (assigned by max entropy)
#The final weight of the representative is the summation of weights for everyone in the cluster
for i = 1:length(unique(T(:,1)))
    aux = find(T(:,1) == i)
    weights(i) = sum(xsol(T(aux,3))) 
    structs(i) = find(xsol == max(xsol(T(aux,3)))) 
end

save weights.mat weights
save structs.mat structs

function [coor,atnam,at_res]=readpdb(fname,reslst,atlst,model,off,chainID) #NEW INPUT FLAG
#-----------------------------------------------------------------
#   df-nov-97
#	read in a pdb data set 
#       uses: readasci.m
#(c) Copyright by David Fushman, U.Maryland
#-----------------------------------------------------------------
if nargin<6, chainID='A'   end
if nargin<5, off=0         end
if nargin<4, model=1       end		#default: first model
if nargin<3, atlst=[]      end      #default: all atoms
if nargin<2, reslst=[]     end		#default: all resid.

pdb=readasci(fname)
#disp(['got the ',fname,' data set, analyzing...'])
nlin=size(pdb,1)
select_at=zeros(nlin,1)
isel=1
nmod=1
termflag=0                             #term flag off
for ii=1:nlin,
  if strcmp('TER',pdb(ii,1:3))|strcmp('END',pdb(ii,1:3))|strcmp('.',pdb(ii,1)), 
     if termflag==0,
        nmod=nmod+1
        termflag=1                     #term flag on
        if nmod>model, break end
     end
  else
     termflag=0                        #term flag off
     if nmod==model                     #read in the structure
       if strcmp('ATOM',pdb(ii,1:4))&(strcmp(pdb(ii,22),chainID)== 1), 
         if ~isempty(reslst),           #select res#
            if ~isempty(find(reslst==str2num(pdb(ii,23:26)))),
               select_at(isel)=ii isel=isel+1
            end
         else  select_at(isel)=ii isel=isel+1
         end
       end
     end
  end
end
sel=select_at(1:isel-1)
nat=length(sel)
if nat==0, error('no atoms found!!! wrong filename or model') end  
#disp([num2str(nat),' atoms read in'])
atnam=pdb(sel,8:26)
coor=zeros(nat,3)
at_res=zeros(nat,2)
at_res(:,1)=str2num(pdb(sel,7:11))
at_res(:,2)=str2num(pdb(sel,23:26))
coor(:,1)=str2num(pdb(sel,31:38))
coor(:,2)=str2num(pdb(sel,39:46))
coor(:,3)=str2num(pdb(sel,47:54))
#coor(:,1:3)=str2num(pdb(sel,31:54))

#select only required atoms 
if ~isempty(atlst),
   rlist=[min(at_res(:,2)):max(at_res(:,2))]
   selat=at_select(atnam,at_res,rlist,atlst,2,off) #NEW INPUTS IN AT_SELECT  
   coor=coor(selat,:)
   at_res=at_res(selat,:)
   atnam=atnam(selat,:)
   disp([num2str(length(selat)),' atoms selected'])
end    
end
#=================================================================
function [coor]=readpdb1(fname,reslst,atlst,model,tag1,tag2) #NEW INPUT FLAG
#-----------------------------------------------------------------
#   df-nov-97
#	read in a pdb data set 
#       uses: readasci.m
#(c) Copyright by David Fushman, U.Maryland
#-----------------------------------------------------------------
if nargin<4, model=1 	end		#default: first model
if nargin<3, atlst=[]  end      #default: all atoms
if nargin<2, reslst=[] end		#default: all resid.

selat = tag1(:,1)
sel = tag2(:,1)+1

pdb=readasci(fname)
#disp(['got the ',fname,' data set, analyzing...'])
nlin=size(pdb,1)
select_at=zeros(nlin,1)
nmod=1

nat=length(sel)
if nat==0, error('no atoms found!!! wrong filename or model') end  
#disp([num2str(nat),' atoms read in'])
coor=zeros(nat,3)
coor(sel-1,1)=str2num(pdb(sel,31:38))
coor(sel-1,2)=str2num(pdb(sel,39:46))
coor(sel-1,3)=str2num(pdb(sel,47:54))
#coor(:,1:3)=str2num(pdb(sel,31:54))

#select only required atoms 
if ~isempty(atlst),
   coor=coor(selat,:)
   #disp([num2str(length(selat)),' atoms selected'])
else
    coor=coor(sel-1,:)
end    

end
#=================================================================
function z=readasci(filename)
#---------------------------------------------------
#  df-feb-08 (took care of EOF problems) df-nov-97
#    read in the an ASCII file using fread
#---------------------------------------------------
fid = fopen(filename,'r')
F = fread(fid)
fclose(fid)
ENDLINE=10
endl=(find(F==ENDLINE))
tot_len=length(F)
if endl(end)~=tot_len, endl=[endltot_len] end  #if no endline at EOF
nlin=length(endl)
maxlen=endl(1)
for ii=2:nlin,
  if maxlen<endl(ii)-endl(ii-1), maxlen=endl(ii)-endl(ii-1)end
end
mtrx=ones(nlin,maxlen)*32
mtrx(1,1:endl(1)-1)=F(1:endl(1)-1)'	#first line
for ii=2:nlin,
  mtrx(ii,1:endl(ii)-endl(ii-1)-1)=F(endl(ii-1)+1:endl(ii)-1)'
end
z=setstr(mtrx)
end
#====================================================
function [coor3s,Rotmat,rmsd,coor2s]=superimpose(coor1,atnam1,at_res1,coor2,atnam2,at_res2,reslst1,reslst2,atsel,coor3)
#---------------------------------------------------
#   df-jun-13  fixed the bug in atnam positions
#   df-dec-09 (takes care of differences in atom order)  
#   df-aug-03
#   given a list of residues and atom selections
#   superimpose coor2 onto coor1, 
#   if coor3 is given, rotate this set as well
#   if not, rotated coor2 will be reported as coor3s
#---------------------------------------------------

#preliminary selection
sel1=at_select(atnam1,at_res1,reslst1,atsel,3)
#sel1=at_select(atnam1,at_res1,reslst1,atsel,2,-1)
sel2=at_select(atnam2,at_res2,reslst2,atsel,3)
#sel2=at_select(atnam2,at_res2,reslst2,atsel,2,-1)

#length(sel1)
#length(sel2)
if length(sel1)~=length(sel2), 
    error('sets of atom coordinates are different!')
end

#matching atom names
sel2m=sel2
for ii=1:length(reslst1),
   ind1=find(at_res1(sel1(:),2)==reslst1(ii))
   ind2=find(at_res2(sel2(:),2)==reslst2(ii))
   if ~isempty(ind1)&~isempty(ind2),
      if length(ind1)~=length(ind2), 
            error(['atoms in residues', num2str(reslst1(ii)),' and ',num2str(reslst2(ii)),' are different!'])
      end 
       #[atnam1(sel1(ind1),:),atnam2(sel2(ind2),:)]
      ind2m=ind2
      for jj=1:length(ind1),
          for kk=1:length(ind2),
             if strcmp(atnam1(sel1(ind1(jj)),6:9),atnam2(sel2(ind2(kk)),6:9))   
                ind2m(jj)=ind2(kk)
                break
             end
          end
      end
      sel2m(ind2)=sel2(ind2m)
      #disp('modif')
      #[atnam1(sel1(ind1),:),atnam2(sel2m(ind2),:)]
   end    
end
        
#sel1=at_select(atnam1,at_res1,reslst1,atsel)
#sel2=at_select(atnam2,at_res2,reslst2,atsel)
#length(sel1)
#length(sel2m)
#[atnam1(sel1,:),atnam2(sel2,:)]
#if length(sel1)~=length(sel2), 
#    error('sets of atom coordinates are different!')
#end
#NEED to coordinate atom names in both sets!!!
if exist('coor3'),
   [coor3s,Rotmat,coor2s,rmsd]=rotfit(coor1(sel1,:),coor2(sel2m,:),coor3)
else
   [coor3s,Rotmat,coor2s,rmsd]=rotfit(coor1(sel1,:),coor2(sel2m,:))
end
#disp(['rmsd= ',num2str(rmsd)])
#===================================================
end

function sel=at_select(atnam,at_res,reslst,atlst,atdim,offs),
#---------------------------------------------------
#  df-aug-03    df-dec-97
#	select atoms according to atnam & reslst
#   standard lists: 'all', 'heavy', 'bb', 'nh'
#   atdim is required if atlist in not a standard one
#---------------------------------------------------
if nargin<6, offs=0 end		#default  
if nargin<5, atdim=2 end		#default
atsel_flag=1
if strcmp(atlst,'heavy')|strcmp(atlst,'HEAVY')
    at='NCOS' atlen=1
elseif strcmp(atlst,'bb')|strcmp(atlst,'BB')
    at='N CAC O ' atlen=2
elseif strcmp(atlst,'nh')|strcmp(atlst,'NH') 
   #not now 
   at='HN' atlen=2
elseif strcmp(atlst,'all')|strcmp(atlst,'ALL')|isempty(atlst), 
   #get all atoms
   #at='HNCOS'atlen=1 
   atsel_flag=0 atlen=1
else
    atlen=atdim
    at=atlst
end

atnampos=[7+offs:7+offs+atlen-1]
nat=size(atnam,1)
indsel=NaN*ones(nat,1)
for ii=1:nat,
  if isempty(reslst),			#take all resid
    if atsel_flag==0,           	#all atoms
        indsel(ii)=ii
    else                        	#select atoms
        if ~isempty(findstr(atnam(ii,atnampos),at)),
          indsel(ii)=ii 
        end
    end
  else					        #res. selection
   if ~isempty(find(at_res(ii,2)==reslst)),
     if atsel_flag==0,			#all atoms
        indsel(ii)=ii 
     else				      #select atoms
        if ~isempty(findstr(atnam(ii,atnampos),at)),
             indsel(ii)=ii   
        end
     end
   end
  end
end
sel=indsel(~isnan(indsel))
end

function [b2rot,rot,b1rot,rmsd]=rotfit(a,b1,b2)
#----------------------------------------------------
# df-aug-03     df-nov-97		df-jul-97	df-jan-93
#    perform a rotation of coord.set <b1>  
#    in order to superimpose it onto <a>
#    and then apply this rotation to coord.set <b2>  
#    as a first step the centers-of-mass (geom.centers)
#    are superimposed
#    it is assumed that <b1> is a subset of <b2>
#
#    for the algorithm see: 
#    A.D.McLachlan, J.Mol.Biol.,1979,v.128,49 (appendix)
#
#    CALL: [b2rot,rot,b1rot]=mdrotfit(a,b1,b2)
#          a -  coordinates (number-of-atoms x 3)
#          b1 - coordinates (number-of-atoms x 3)
#	   b2 - coordinates (number-of-atoms x 3)
#    OUTPUT: b2rot, b1rot being the result of a rotation of b2
#	     rot is a rotation matrix, superimposing 
#	     <b1> onto <a>
#    DOESN'T work for number-of-atoms = 1  !!!
#    Modification of nov97: removed row2coor, mtrx2row
#----------------------------------------------------
 sets=2
 if nargin<3, b2 = [] end		#only one set to rotate
 if isempty(b2), sets=1 end    #only one set to rotate
 nat=size(a,1) 
 natb1=size(b1,1) 
 if nat~=natb1, error('different number of atoms in the structures!') end
 if nat==1, error('DOESN"T work for number-of-atoms = 1 !!!') end
 if sets==2, natb2=size(b2,1) end 
 cma=mean(a)
 cmb1=mean(b1)
 # shift to put c-o-mass in origin
 as=a-ones(nat,1)*cma
 b1s=b1-ones(nat,1)*cmb1
 if sets==2, b2s=b2-ones(natb2,1)*cmb1 end
 detind=1
 u=zeros(3,3)
 omega=zeros(6,6)
 u=as'*b1s
 if det(u)<0, disp('Determinant negative!!!') detind=0 end   # negative-det(u)-records
 omega(1:3,4:6)=u				#build omega-matrix
 omega(4:6,1:3)=u'
 [V,D]=eig(omega)              	#eigenvalues of omega
 rot=zeros(3,3) 
 H1=zeros(3,3)
 K1=zeros(3,3)
 r1=zeros(3,3)
 for k=1:6,
   K1=ones(3,1)*V(1:3,k)'
   H1=V(4:6,k)*ones(1,3)
   r1=K1.*H1
   rot=rot+sign(D(k,k))*r1 
 end
 b1rot=b1s*rot				#(rot'*b1s')'	
 #calculate RMSDs
 diff=as-b1rot
 rmsd=sqrt(mean(sum((diff.^2)')'))
 #shift <b>'s to COM of <a>
 b1rot=b1rot+ones(nat,1)*cma		# shift back to c-o-mass of <a>
 if sets==2, 
   b2rot=b2s*rot			#(rot'*b2s')'
   b2rot=b2rot+ones(natb2,1)*cma 	# shift back to c-o-mass of <a>
 else
   b2rot=b1rot
 end
end
