//////////////////////////////////////////////////////////////////////////
//
// gHBfix program, version 1.3
// original version of this file was attached as Supporting Information to the article
// Kuhrova, P., Mlynsky, V., et al. IMPROVING THE PERFORMANCE OF THE RNA
// AMBER FORCE FIELD BY TUNING THE HYDROGEN-BONDING INTERACTIONS,
// J. Chem. Theory Comput. 2019, 15, 5, 3288–3305, (https://doi.org/10.1021/acs.jctc.8b00955)
// gHBfix.cpp is the only source code for gHBfix program.
// Copyright (C) 2022  Pavel Banas, Regional Centre of Advanced Technologies and Materials,
//                     Czech Advanced Technology and Research Institute (CATRIN),
//                     Palacky University Olomouc, Czech Republic
//                     pavel(dot)banas(et)upol(dot)cz
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.
//
//////////////////////////////////////////////////////////////////////////

#include <cstring>
#include <cstdio>
#include <vector>
#include <cmath>

using std::vector;

class nt
{
public:
  nt();
  ~nt();
  void Init(char base, int last, bool ter5=false, bool ter3=false);
  bool O2O3test(int O2, int O3);

public:
 int         m_atoms;
 vector<int> m_don_NH;
 vector<int> m_don_OH;
 vector<int> m_acc_Nbase;
 vector<int> m_acc_Obase;
 vector<int> m_acc_OH;
 vector<int> m_acc_nbO;
 vector<int> m_acc_bO;
 vector<int> m_acc_O4;
};

nt::nt()
{
  m_atoms=0;
}

nt::~nt()
{
}

void nt::Init(char base, int last, bool ter5, bool ter3)
{
  //register OH donors, keep O2' always first to check O2'...O3' pair within the nt (then O3' if present and after that O5')
  int h2 = last+18;
  if(ter5==true) h2-=2;
  switch(base)
  {
    case 'g':
    case 'G':
      h2+=15;
      break;
    case 'a':
    case 'A':
      h2+=14;
      break;
    case 'u':
    case 'U':
      h2+=11;
      break;
    case 'c':
    case 'C':
      h2+=12;
      break;
    case 'p':
    case 'P':
      h2+=13;
      break;
    case 'm':
    case 'M':
      h2+=15;
      break;
    default:
      printf("Error: inccorect base name in nt::Init routine\n");
  }
  m_don_OH.push_back(h2); // HO2'
  m_acc_OH.push_back(h2-1); // O2'
  if(ter3==true)
  {
    m_don_OH.push_back(h2+2);  // H3T
    m_acc_OH.push_back(h2+1);  // O3' (in case of O3'-H3T)
  }
  if(ter5==true)
  {
    m_don_OH.push_back(last+1); // H5T
    m_acc_OH.push_back(last+2); // O5' (in case of O5'-H5T)
  }

  //register phosphate acceptors
  if(ter5==false)
  {
    m_acc_nbO.push_back(last+2); // O1P
    m_acc_nbO.push_back(last+3); // O2P
    m_acc_bO.push_back(last+4); // O5'
  }
  if(ter3==false) m_acc_bO.push_back(h2+1); // O3'

  //register O4'
  if(ter5==true) m_acc_O4.push_back(last+8);
  else m_acc_O4.push_back(last+10);

  //register base donors and acceptors
  int h1=last+12;
  if(ter5==true) h1-=2;
  switch(base)
  {
    case 'g':
    case 'G':
      m_don_NH.push_back(h1+9); // H1
      m_don_NH.push_back(h1+12); // H21
      m_don_NH.push_back(h1+13); // H22
      m_acc_Nbase.push_back(h1+4); // N7
      m_acc_Obase.push_back(h1+7); // O6
      m_acc_Nbase.push_back(h1+14); // N3
      break;
    case 'a':
    case 'A':
      m_don_NH.push_back(h1+8); // H61
      m_don_NH.push_back(h1+9); // H62
      m_acc_Nbase.push_back(h1+4); // N7
      m_acc_Nbase.push_back(h1+10); // N1
      m_acc_Nbase.push_back(h1+13); // N3
      break;
    case 'c':
    case 'C':
      m_don_NH.push_back(h1+8); // H41
      m_don_NH.push_back(h1+9); // H42
      m_acc_Nbase.push_back(h1+10); // N3
      m_acc_Obase.push_back(h1+12); // O2
      break;
    case 'u':
    case 'U':
      m_don_NH.push_back(h1+9); // H3
      m_acc_Obase.push_back(h1+7); // O2
      m_acc_Obase.push_back(h1+11); // O4
      break;
    case 'p':
    case 'P':
      m_don_NH.push_back(h1+8); // H41
      m_don_NH.push_back(h1+9); // H42
      m_don_NH.push_back(h1+11); // H3
      m_acc_Obase.push_back(h1+13); // O2
      break;
    case 'm':
    case 'M':
      m_don_NH.push_back(h1+8); // H61
      m_don_NH.push_back(h1+9); // H62
      m_don_NH.push_back(h1+11); // H1
      m_acc_Nbase.push_back(h1+4); // N7
      m_acc_Nbase.push_back(h1+14); // N3
      break;
    default:
      printf("Error: inccorect base name in nt::Init routine\n");
  }

  //set m_atoms
  switch(base)
  {
    case 'g':
    case 'G':
      m_atoms=34;
      break;
    case 'a':
    case 'A':
      m_atoms=33;
      break;
    case 'c':
    case 'C':
      m_atoms=31;
      break;
    case 'u':
    case 'U':
      m_atoms=30;
      break;
    case 'p':
    case 'P':
      m_atoms=32;
      break;
    case 'm':
    case 'M':
      m_atoms=34;
      break;
    default:
      printf("Error: inccorect base name in nt::Init routine\n");
  }
  if(ter5==true) m_atoms-=2;
  if(ter3==true) m_atoms+=1;

  return;
}

bool nt::O2O3test(int HO2, int O3)
{
  int i;
  for(i=0;i<m_don_OH.size();i++) if((m_don_OH[i]==HO2)&&(m_acc_OH[i]==O3)) return true; // check -OH self interaction
  if(((m_acc_bO.size()+m_acc_nbO.size())==1)||((m_acc_bO.size()+m_acc_nbO.size())==4)) // O3' is bO acceptor
  {
    if((m_don_OH[0]==HO2)&&(m_acc_bO[m_acc_bO.size()-1]==O3)) return true;
  }
  else // O3' is OH occeptor and its H3T is OH donor, so we have to check both O3'-H...O2' and O2'-H...O3'
  {
    if((m_don_OH[0]==HO2)&&(m_acc_OH[1]==O3)) return true;
    if((m_don_OH[1]==HO2)&&(m_acc_OH[0]==O3)) return true;
  }
  return false;
}

void PrintHBfix(FILE* output, int i, int j, double eps,  const char* c, int* n, int format, vector<int>* p_donors, vector<int>* p_acceptors)
{
  // format 0: AMBER, 1: GROMACS(plumed) (currently disabled in this code), 2: cpptraj, 3: GROMACS(plumed) gHBfix format
  if(fabs(eps)>0.0000001)
  {
    (*n)++;
    if(format==0)
    {
      fprintf(output,"&rst iat=%d, %d,\n",i,j);
      fprintf(output,"iresid=0,nstep1=0,nstep2=0,\nirstyp=0,ifvari=0,ninc=0,imult=0,ir6=0,ifntyp=0,\n");
      fprintf(output,"r1=0.0, r2=2.0, r3=2.0, r4=3.0, rk2=0.0000, rk3=%.10lf,\n/\n",-5.0*eps);
      fprintf(output,"&rst iat=%d, %d,\n",i,j);
      fprintf(output,"iresid=0,nstep1=0,nstep2=0,\nirstyp=0,ifvari=0,ninc=0,imult=0,ir6=0,ifntyp=0,\n");
      fprintf(output,"r1=0.0, r2=2.0, r3=2.0, r4=2.8, rk2=0.0000, rk3=%.10lf,\n/\n",6.25*eps);
    }
    else if(format==1)
    {
      fprintf(output,"%s%d: DISTANCE ATOMS=%d,%d\n",c,*n,i,j);
      fprintf(output,"f_%s%da: MATHEVAL ARG=%s%d FUNC=((-%.10lf)*(0.2-x)^2)*step(x-0.2)*step(0.3-x)+((-%.10lf)*(0.1)*(2*x-0.5))*step(x-0.3) PERIODIC=NO\n",c,*n,c,*n,5.0*eps*418.4,5.0*eps*418.4);
      fprintf(output,"f_%s%db: MATHEVAL ARG=%s%d FUNC=((%.10lf)*(0.2-x)^2)*step(x-0.2)*step(0.28-x)+((%.10lf)*(0.08)*(2*x-0.48))*step(x-0.28) PERIODIC=NO\n",c,*n,c,*n,6.25*eps*418.4,6.25*eps*418.4);
    }
    else if(format==2)
    {
      fprintf(output,"distance %s%d @%d @%d\n",c,*n,i,j);
    }
    else if(format==3) // no printing, just save the indexes and print them later together
    {
      p_donors->push_back(i);
      p_acceptors->push_back(j);
    }
    else printf("Error: incorrect output format in PrintHBfix routine.\n");
  }
  return;
}

int main(int argc, char* argv[])
{
  if(argc<3)
  {
    printf(" usage: %s <sequence> <output_file> [ -lambda <lambda> ] [ -hb-<don1>-<acc1> <eta1_value> -hb-<don2>-<acc2> <eta2_value> ] [ -plumed ] [ -cpptraj ]\n\nwrite sequence in the following form GCGCGC,GCGCGC (RNA strands separated by comma without any spaces, use P for protonated Cyt and M for protonated Ade)\ndonors: NH, OH\nacceptors: nbO, bO, OH (oxygens of OH groups), O4, O (base oxygens), N (base nitrogens)\nsuggested gHBfix setting: -hb-NH-N 1.0 -hb-OH-nbO -0.5 (-hb-OH-bO -0.5)\n-plumed writes gHBfix in plumed gHBfix format\n-cpptraj writes input for cpptraj rather than gHBfix input file\n",argv[0]);
//    printf("number of inputs: %d\n",argc);
//    for(int ii=0; ii<argc;ii++) printf("arg %d: %s\n",argc,argv[ii]);
    return 0;
  }

  //parse command line
  int i,j;
  int last=0;
  nt* p_nt;
  vector<nt*> seq;
  bool ter5=true,ter3=false;
  for(i=0;i<strlen(argv[1]);i++)
  {
    if((argv[1][i]=='G')||(argv[1][i]=='g')||(argv[1][i]=='A')||(argv[1][i]=='a')||(argv[1][i]=='C')||(argv[1][i]=='c')||(argv[1][i]=='U')||(argv[1][i]=='u')||(argv[1][i]==',')||(argv[1][i]=='P')||(argv[1][i]=='p')||(argv[1][i]=='M')||(argv[1][i]=='m'))
    {
      //set variables for terminal nts
      if(i!=0)
      {
        if(argv[1][i-1]==',') ter5=true;
        else ter5=false;
      }
      if(i==strlen(argv[1])-1) ter3=true;
      else if(argv[1][i+1]==',') ter3=true;
      else ter3=false;

      //register nt
      if(argv[1][i]!=',')
      {
        p_nt = new nt;
        p_nt->Init(argv[1][i],last,ter5,ter3);
        last+=p_nt->m_atoms;
        seq.push_back(p_nt);
      }
    }
    else
    {
      printf("Error: incorrect character in the sequence %s",argv[1]);
      printf(" usage: %s <sequence> <output_file> [ -lambda <lambda> ] [ -hb-<don1>-<acc1> <eta1_value> -hb-<don2>-<acc2> <eta2_value> ] [ -plumed ] [ -cpptraj ]\n\nwrite sequence in the following form GCGCGC,GCGCGC (RNA strands separated by comma without any spaces, use P for protonated Cyt and M for protonated Ade)\ndonors: NH, OH\nacceptors: nbO, bO, OH (oxygens of OH groups), O4, O (base oxygens), N (base nitrogens)\nsuggested gHBfix setting: -hb-NH-N 1.0 -hb-OH-nbO -0.5 (-hb-OH-bO -0.5)\n-plumed writes gHBfix in plumed gHBfix format\n-cpptraj writes input for cpptraj rather than gHBfix input file\n",argv[0]);
      return -1;
    }
  }
  double eta[2][6];
  char* pch;
  int don,acc; // don(0:NH, 1:OH), acc(0:nbO, 1:bO, 2:OH, 3:O4, 4:O, 5:N)
  for(i=0;i<2;i++)
  {
    for(j=0;j<6;j++) eta[i][j]=0.0;
  }
  double lambda=1.0;
  bool cpptraj=false;
  bool plumed=false;
  for(i=3;i<argc;i++)
  {
    if((strcmp(argv[i],"-lambda")==0)&&(i+1<argc))
    {
      i++;
      if(sscanf(argv[i],"%lf",&lambda)!=1)
      {
        printf("Error: cannot read %s value\n",argv[i-1]);
        return -1;
      }
    }
    else if((strncmp(argv[i],"-hb-",4)==0)&&(i+1<argc))
    {
      pch=strtok(&(argv[i][3]),"-");
      if(strcmp(pch,"NH")==0) don=0;
      else if(strcmp(pch,"OH")==0) don=1;
      else
      {
        printf("Error: incorrect donor in %s\n",argv[i]);
        return -1;
      }
      pch=strtok(NULL,"-");
      if(strcmp(pch,"nbO")==0) acc=0;
      else if(strcmp(pch,"bO")==0) acc=1;
      else if(strcmp(pch,"OH")==0) acc=2;
      else if(strcmp(pch,"O4")==0) acc=3;
      else if(strcmp(pch,"O")==0) acc=4;
      else if(strcmp(pch,"N")==0) acc=5;
      else
      {
        printf("Error: incorrect acceptor in %s\n",argv[i]);
        return -1;
      }
      i++;
      if(sscanf(argv[i],"%lf",&(eta[don][acc]))!=1)
      {
        printf("Error: cannot read %s value\n",argv[i-1]);
        return -1;
      }
    }
    else if(strcmp(argv[i],"-plumed")==0) plumed=true;
    else if(strcmp(argv[i],"-cpptraj")==0) cpptraj=true;
    else
    {
      printf("Error: unexpected input in parsing command line\n");
      printf(" usage: %s <sequence> <output_file> [ -lambda <lambda> ] [ -hb-<don1>-<acc1> <eta1_value> -hb-<don2>-<acc2> <eta2_value> ] [ -plumed ] [ -cpptraj ]\n\nwrite sequence in the following form GCGCGC,GCGCGC (RNA strands separated by comma without any spaces, use P for protonated Cyt and M for protonated Ade)\ndonors: NH, OH\nacceptors: nbO, bO, OH (oxygens of OH groups), O4, O (base oxygens), N (base nitrogens)\nsuggested gHBfix setting: -hb-NH-N 1.0 -hb-OH-nbO -0.5 (-hb-OH-bO -0.5)\n-plumed writes gHBfix in plumed gHBfix format\n-cpptraj writes input for cpptraj rather than gHBfix input file\n",argv[0]);
      return -1;
    }
  }
  for(i=0;i<2;i++)
  {
    for(j=0;j<6;j++) eta[i][j]*=lambda;
  }

  // print all HBfix restraints
  FILE* output;
  int output_format=0; // 0: AMBER, 1: GROMACS(plumed), 2: cpptraj, 3: GROMACS(plumed) gHBfix format
  if(plumed==true) output_format=3;
  if(cpptraj==true) output_format=2;
  if((output = fopen(argv[2],"w"))==NULL)
  {
    printf("Error: cannot open output file %s\n",argv[2]);
    return -1;
  }
  if(output_format==1) fprintf(output,"#WARNING: AMBER atom ordering is assumed, for gromacs this has to be checked!!\n");
  vector<int> donors;
  vector<int> acceptors;
  int k,l;
  int na,nb,nc,nd,ne,nf,ng,nh,ni,nj,nk,nl;
  na=nb=nc=nd=ne=nf=ng=nh=ni=nj=nk=nl=0;
  for(i=0;i<seq.size();i++)
  {
    for(j=0;j<seq.size();j++)
    {
      // NH donors
      for(k=0;k<seq[i]->m_don_NH.size();k++)
      {
        // nbO acceptors
        for(l=0;l<seq[j]->m_acc_nbO.size();l++)
          PrintHBfix(output,seq[i]->m_don_NH[k],seq[j]->m_acc_nbO[l],eta[0][0],"a",&na,output_format,&donors,&acceptors);
        // bO acceptors
        for(l=0;l<seq[j]->m_acc_bO.size();l++)
          PrintHBfix(output,seq[i]->m_don_NH[k],seq[j]->m_acc_bO[l],eta[0][1],"b",&nb,output_format,&donors,&acceptors);
        // OH acceptors
        for(l=0;l<seq[j]->m_acc_OH.size();l++)
          PrintHBfix(output,seq[i]->m_don_NH[k],seq[j]->m_acc_OH[l],eta[0][2],"c",&nc,output_format,&donors,&acceptors);
        // O4 acceptors
        for(l=0;l<seq[j]->m_acc_O4.size();l++)
          PrintHBfix(output,seq[i]->m_don_NH[k],seq[j]->m_acc_O4[l],eta[0][3],"d",&nd,output_format,&donors,&acceptors);
        if(i==j) continue; //skip intrabase NH...N and NH...O interactions
        // Obase acceptors
        for(l=0;l<seq[j]->m_acc_Obase.size();l++)
          PrintHBfix(output,seq[i]->m_don_NH[k],seq[j]->m_acc_Obase[l],eta[0][4],"e",&ne,output_format,&donors,&acceptors);
        // Nbase acceptors
        for(l=0;l<seq[j]->m_acc_Nbase.size();l++)
          PrintHBfix(output,seq[i]->m_don_NH[k],seq[j]->m_acc_Nbase[l],eta[0][5],"f",&nf,output_format,&donors,&acceptors);
      }
      // OH donors
      for(k=0;k<seq[i]->m_don_OH.size();k++)
      {
        // nbO acceptors
        for(l=0;l<seq[j]->m_acc_nbO.size();l++)
          PrintHBfix(output,seq[i]->m_don_OH[k],seq[j]->m_acc_nbO[l],eta[1][0],"g",&ng,output_format,&donors,&acceptors);
        // bO acceptors
        for(l=0;l<seq[j]->m_acc_bO.size();l++)
        {
          if((i==j)&&(seq[i]->O2O3test(seq[i]->m_don_OH[k],seq[j]->m_acc_bO[l])==true)) continue; //skip O2'...O3' interaction within the nucleotide
          PrintHBfix(output,seq[i]->m_don_OH[k],seq[j]->m_acc_bO[l],eta[1][1],"h",&nh,output_format,&donors,&acceptors);
        }
        // OH acceptors
        for(l=0;l<seq[j]->m_acc_OH.size();l++)
        {
          if((i==j)&&(seq[i]->O2O3test(seq[i]->m_don_OH[k],seq[j]->m_acc_OH[l])==true)) continue; //skip O2'...O3' interaction within the nucleotide
          PrintHBfix(output,seq[i]->m_don_OH[k],seq[j]->m_acc_OH[l],eta[1][2],"i",&ni,output_format,&donors,&acceptors);
        }
        // O4 acceptors
        for(l=0;l<seq[j]->m_acc_O4.size();l++)
          PrintHBfix(output,seq[i]->m_don_OH[k],seq[j]->m_acc_O4[l],eta[1][3],"j",&nj,output_format,&donors,&acceptors);
        // Obase acceptors
        for(l=0;l<seq[j]->m_acc_Obase.size();l++)
          PrintHBfix(output,seq[i]->m_don_OH[k],seq[j]->m_acc_Obase[l],eta[1][4],"k",&nk,output_format,&donors,&acceptors);
        // Nbase acceptors
        for(l=0;l<seq[j]->m_acc_Nbase.size();l++)
          PrintHBfix(output,seq[i]->m_don_OH[k],seq[j]->m_acc_Nbase[l],eta[1][5],"l",&nl,output_format,&donors,&acceptors);
      }
    }
  }
  if(output_format==2)
  {
    fprintf(output,"run\n");
    fprintf(output,"writedata don_NH-acc_nbO.dat ");
    for(i=1;i<=na;i++) fprintf(output,"a%d ",i);
    fprintf(output,"\nwritedata don_NH-acc_bO.dat ");
    for(i=1;i<=nb;i++) fprintf(output,"b%d ",i);
    fprintf(output,"\nwritedata don_NH-acc_OH.dat ");
    for(i=1;i<=nc;i++) fprintf(output,"c%d ",i);
    fprintf(output,"\nwritedata don_NH-acc_O4.dat ");
    for(i=1;i<=nd;i++) fprintf(output,"d%d ",i);
    fprintf(output,"\nwritedata don_NH-acc_Obase.dat ");
    for(i=1;i<=ne;i++) fprintf(output,"e%d ",i);
    fprintf(output,"\nwritedata don_NH-acc_Nbase.dat ");
    for(i=1;i<=nf;i++) fprintf(output,"f%d ",i);
    fprintf(output,"\nwritedata don_OH-acc_nbO.dat ");
    for(i=1;i<=ng;i++) fprintf(output,"g%d ",i);
    fprintf(output,"\nwritedata don_OH-acc_bO.dat ");
    for(i=1;i<=nh;i++) fprintf(output,"h%d ",i);
    fprintf(output,"\nwritedata don_OH-acc_OH.dat ");
    for(i=1;i<=ni;i++) fprintf(output,"i%d ",i);
    fprintf(output,"\nwritedata don_OH-acc_O4.dat ");
    for(i=1;i<=nj;i++) fprintf(output,"j%d ",i);
    fprintf(output,"\nwritedata don_OH-acc_Obase.dat ");
    for(i=1;i<=nk;i++) fprintf(output,"k%d ",i);
    fprintf(output,"\nwritedata don_OH-acc_Nbase.dat ");
    for(i=1;i<=nl;i++) fprintf(output,"l%d ",i);
    fprintf(output,"\n");
  }
  else if(output_format==1)
  {
    if(na>0)
    {
      fprintf(output,"BIASVALUE ARG=");
      for(i=1;i<na;i++) fprintf(output,"f_a%da,f_a%db,",i,i);
      fprintf(output,"f_a%da,f_a%db",na,na);
      fprintf(output," LABEL=b_a\n");
    }
    if(nb>0)
    {
      fprintf(output,"BIASVALUE ARG=");
      for(i=1;i<nb;i++) fprintf(output,"f_b%da,f_b%db,",i,i);
      fprintf(output,"f_b%da,f_b%db",nb,nb);
      fprintf(output," LABEL=b_b\n");
    }
    if(nc>0)
    {
      fprintf(output,"BIASVALUE ARG=");
      for(i=1;i<nc;i++) fprintf(output,"f_c%da,f_c%db,",i,i);
      fprintf(output,"f_c%da,f_c%db",nc,nc);
      fprintf(output," LABEL=b_c\n");
    }
    if(nd>0)
    {
      fprintf(output,"BIASVALUE ARG=");
      for(i=1;i<nd;i++) fprintf(output,"f_d%da,f_d%db,",i,i);
      fprintf(output,"f_d%da,f_d%db",nd,nd);
      fprintf(output," LABEL=b_d\n");
    }
    if(ne>0)
    {
      fprintf(output,"BIASVALUE ARG=");
      for(i=1;i<ne;i++) fprintf(output,"f_e%da,f_e%db,",i,i);
      fprintf(output,"f_e%da,f_e%db",ne,ne);
      fprintf(output," LABEL=b_e\n");
    }
    if(nf>0)
    {
      fprintf(output,"BIASVALUE ARG=");
      for(i=1;i<nf;i++) fprintf(output,"f_f%da,f_f%db,",i,i);
      fprintf(output,"f_f%da,f_f%db",nf,nf);
      fprintf(output," LABEL=b_f\n");
    }
    if(ng>0)
    {
      fprintf(output,"BIASVALUE ARG=");
      for(i=1;i<ng;i++) fprintf(output,"f_g%da,f_g%db,",i,i);
      fprintf(output,"f_g%da,f_g%db",ng,ng);
      fprintf(output," LABEL=b_g\n");
    }
    if(nh>0)
    {
      fprintf(output,"BIASVALUE ARG=");
      for(i=1;i<nh;i++) fprintf(output,"f_h%da,f_h%db,",i,i);
      fprintf(output,"f_h%da,f_h%db",nh,nh);
      fprintf(output," LABEL=b_h\n");
    }
    if(ni>0)
    {
      fprintf(output,"BIASVALUE ARG=");
      for(i=1;i<ni;i++) fprintf(output,"f_i%da,f_i%db,",i,i);
      fprintf(output,"f_i%da,f_i%db",ni,ni);
      fprintf(output," LABEL=b_i\n");
    }
    if(nj>0)
    {
      fprintf(output,"BIASVALUE ARG=");
      for(i=1;i<nj;i++) fprintf(output,"f_j%da,f_j%db,",i,i);
      fprintf(output,"f_j%da,f_j%db",nj,nj);
      fprintf(output," LABEL=b_j\n");
    }
    if(nk>0)
    {
      fprintf(output,"BIASVALUE ARG=");
      for(i=1;i<nk;i++) fprintf(output,"f_k%da,f_k%db,",i,i);
      fprintf(output,"f_k%da,f_k%db",nk,nk);
      fprintf(output," LABEL=b_k\n");
    }
    if(nl>0)
    {
      fprintf(output,"BIASVALUE ARG=");
      for(i=1;i<nl;i++) fprintf(output,"f_l%da,f_l%db,",i,i);
      fprintf(output,"f_l%da,f_l%db",nl,nl);
      fprintf(output," LABEL=b_l\n");
    }
  }
  else if(output_format==3)
  {
    fprintf(output,"group_1: GHBFIX PAIR GROUPA=");
    for(i=0;i<donors.size()-1;i++) fprintf(output,"%d,",donors[i]);
    fprintf(output,"%d GROUPB=",donors[donors.size()-1]);
    for(i=0;i<acceptors.size()-1;i++) fprintf(output,"%d,",acceptors[i]);
    fprintf(output,"%d D_0=0.2 D_MAX=0.3 C=0.8 TYPES=typesTable.dat PARAMS=scalingParameters.dat\n",acceptors[acceptors.size()-1]);
    // now create typesTable.dat and scalingParameters.dat
    int index[8];
    for(i=0;i<8;i++) index[i]=0;
    for(i=0;i<2;i++) for(j=0;j<6;j++) if(fabs(eta[i][j])>0.0000001) {index[i]=index[j+2]=1;}
    int num=1;
    for(i=0;i<8;i++) if(index[i]==1) index[i]=num++;
//    for(i=0;i<8;i++) printf("%d\n",index[i]);
    int n_atoms=0;
    for(i=0;i<seq.size();i++) n_atoms+=seq[i]->m_atoms;
    int* atom_index;
    atom_index = new int[n_atoms];
    for(i=0;i<n_atoms;i++) atom_index[i]=0;
    for(i=0;i<seq.size();i++)
    {
      for(j=0;j<seq[i]->m_don_NH.size();j++) atom_index[seq[i]->m_don_NH[j]-1]=index[0];
      for(j=0;j<seq[i]->m_don_OH.size();j++) atom_index[seq[i]->m_don_OH[j]-1]=index[1];
      for(j=0;j<seq[i]->m_acc_nbO.size();j++) atom_index[seq[i]->m_acc_nbO[j]-1]=index[2];
      for(j=0;j<seq[i]->m_acc_bO.size();j++) atom_index[seq[i]->m_acc_bO[j]-1]=index[3];
      for(j=0;j<seq[i]->m_acc_OH.size();j++) atom_index[seq[i]->m_acc_OH[j]-1]=index[4];
      for(j=0;j<seq[i]->m_acc_O4.size();j++) atom_index[seq[i]->m_acc_O4[j]-1]=index[5];
      for(j=0;j<seq[i]->m_acc_Obase.size();j++) atom_index[seq[i]->m_acc_Obase[j]-1]=index[6];
      for(j=0;j<seq[i]->m_acc_Nbase.size();j++) atom_index[seq[i]->m_acc_Nbase[j]-1]=index[7];
    }
    FILE* index_file;
    index_file=fopen("typesTable.dat","w");
    for(i=0;i<n_atoms;i++) fprintf(index_file,"%d\n",atom_index[i]);
    fclose(index_file);
    FILE* param_file;
    param_file=fopen("scalingParameters.dat","w");
    fprintf(param_file,"#! FIELDS itype jtype eta\n");
    for(i=0;i<2;i++) for(j=0;j<6;j++) if(fabs(eta[i][j])>0.0000001) fprintf(param_file,"%d %d %lf\n%d %d %lf\n",index[i],index[j+2],eta[i][j],index[j+2],index[i],eta[i][j]);
    fclose(param_file);
  }
  fclose(output);

  //delete all nts
  for(i=0;i<seq.size();i++) delete seq[i];
  seq.clear();
  return 0;
}
