/*************************************************************************

   Program:    findcore_Apr16
   File:       findcore_Apr16.c
   
   Version:    V1.8
   Date:       05.11.25
   Function:   Find core from multiple structures given the CORA alignment
               file as a staring point
   
   Copyright:  (c) Prof. Andrew C. R. Martin, UCL 1996-2025
   Author:     Prof. Andrew C. R. Martin
   Modified:   Gabby Marsden (nee Reeves)
   Address:    Biomolecular Structure & Modelling Unit,
               Department of Biochemistry & Molecular Biology,
               University College,
               Gower Street,
               London.
               WC1E 6BT.
   EMail:      andrew@bioinf.org.uk
               
**************************************************************************

   This program is released under the GPL V3

**************************************************************************

   Description: 
   ============
   Finds the core of proteins multiply aligned using the CORA program.

**************************************************************************

   Usage:
   ======

**************************************************************************

   Revision History:
   =================
   V1.0  14.11.96 Original
   V1.1  06.12.96 Added -i option for excluding >3A sections
   V1.2  23.01.97 Added -n option to include non-E/H regions which match
   V1.3  13.03.97 Fixed a problem where merging zones could lead to the
                  merged zones having different numbers of residues.
                  e.g. 24-43:24-43 + 41-51:42-52 ==> 24-51:24-52
                  Fixed this by extending zones only if a residue
                  wasn't already in a zone.
   V1.4...V1.7  16.04.02
   V1.8  05.11.25 Updated for new biolib

*************************************************************************/
/* Includes
*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "bioplib/SysDefs.h"
#include "bioplib/MathType.h"
#include "bioplib/pdb.h"
#include "bioplib/general.h"
#include "bioplib/macros.h"
#include "bioplib/fit.h"
#include "bioplib/fsscanf.h"

/************************************************************************/
/* Defines and macros
*/
#define MAXITER 1000
#define MAXBUFF 160
#define MAXCHAR 10
#define DEFAULT_CUT ((REAL)3.0)
#define MAXMALNPNO 50
#define COMMENT

#define TEST(obj)  if( obj == NULL ) printf("no memory for obj !\n")
#define ALLOC(x,y) do { (x) = (y *)malloc(sizeof(y)); TEST(x); } while(0)

typedef struct _zone
{
   struct _zone *next,
                *prev;
   int start[MAXMALNPNO], 
       end[MAXMALNPNO];
}  ZONE;

/*  CORA - Line data for multiple alignment files                       */
typedef struct
{  /* data for each protein in the alignment                            */
   char acid,
        secstruct,
        pdb[6];  
}  Protdata;

typedef struct
{
   int      alnpos,
            conpos,
            proaln,
            colx,
            coly,
            score;
   char     consecstruc;
   /* Pointer defined to point to the first address in Protdata         */
   Protdata *protdata_ptr;   
}  Malndata;

/*  general data for multiple alignment files                           */
typedef struct
{
   int       procnt,
             length; 		                    
   char      title[50],
             proname[MAXMALNPNO][MAXCHAR];
   /* Pointer defined to point to the first address in Malndata         */
   Malndata  *malndata_ptr;      
}  Malign;

/************************************************************************/
/* Globals
*/
BOOL gVerbose      = FALSE,
     gInitialCut   = FALSE,
     gDoRandomCoil = FALSE,
     gDoOutput     = FALSE;

/************************************************************************/
/* Prototypes
*/
int main(int argc, char **argv);
BOOL ParseCmdLine(int argc, char **argv, char *corafile, REAL *dcut);
Malign *ReadCORA(FILE *fp);
ZONE *calcZone(Malign *maln_ptr);
BOOL DefineCore(PDB **pdb, ZONE *zones, Malign *maln_ptr, REAL dcut);
void UpdateBValues(PDB **idx[MAXMALNPNO], int *natoms, int nstruc,
                   ZONE *zones, REAL cutsq);
void SetBValByZone(PDB *pdb, ZONE *zones, int protNum);
BOOL FitCaPDBBFlag(PDB *ref_pdb, PDB *fit_pdb, REAL rm[3][3]);
int CountCore(PDB *pdb);
PDB *DupeCAByBVal(PDB *pdb);
Malign *new_Malign(void);
void clear_Malign(Malign *m);
void WriteTextOutput(ZONE *zones, int *numProts);
void Usage(void);
BOOL DoCut(PDB **idx[MAXMALNPNO], int *natoms, int nstruc,
           ZONE *zones, REAL cutsq);
ZONE *MergeZones(ZONE *zones, int numProts);
BOOL SubsetZone(ZONE *z, ZONE *zones);

/************************************************************************/
int main(int argc, char **argv)
{
   char   corafile[MAXBUFF];
   FILE   *corafp;
   FILE   *pdbfp;
   Malign *maln_ptr;
   int    numPdb;
   int    natoms[MAXMALNPNO];
   int    numProts;
   ZONE   *zones;
   PDB    *pdb[MAXMALNPNO];
   REAL   dcut = DEFAULT_CUT;
   
   if(ParseCmdLine(argc, argv, corafile, &dcut))
   {
      /* Open files                                                     */
      if((corafp=fopen(corafile,"r"))==NULL)
      {
         fprintf(stderr,"Unable to open %s for reading\n",corafile);
         return(1);
      }
      
      if((maln_ptr = ReadCORA(corafp)) == NULL)
      {
         fprintf(stderr,"No Zones read from %s\n",corafile);
         return(1);
      }
      fclose(corafp);
      
      /* calculate inital zones                                         */
      zones = calcZone(maln_ptr);
      numProts = maln_ptr->procnt;
      
      /* open and read the pdbfiles by taking the names from the cora
         file
      */
      for(numPdb = 0; numPdb<maln_ptr->procnt; numPdb++)
      {
         if((pdbfp=fopen(maln_ptr->proname[numPdb],"r"))==NULL)
         {
            fprintf(stderr,"Unable to open %s for reading\n",
                    maln_ptr->proname[numPdb]);
            return(1);
         }
         
         /* Read them and get the data                                  */
         if((pdb[numPdb] = blReadPDB(pdbfp, &natoms[numPdb])) == NULL)
         {
            fprintf(stderr,"No atoms read from PDB file: %s\n",
                    maln_ptr->proname[numPdb]);
            return(1);
         }
         
         fclose(pdbfp);
      }

      /* Print the current zones if required                            */
      if(gVerbose)
      {
         printf("SSAP Zones:\n");
         WriteTextOutput(zones, &numProts);
      }
      
      /* Now call the routine to do the core definition                 */
      DefineCore(pdb, zones, maln_ptr, dcut);
      
      if(gVerbose)
      {
         printf("\nCore before zone merging:\n");
         WriteTextOutput(zones, &numProts);
      }
      
      /* Now remove any zones which are subsets of other zones and merge
         overlapping zones
      */
      zones = MergeZones(zones, numProts);
      
      /* Finally write the output file which lists residues in the 
         structural core and optionally write PDB files with the cores
         flagged
      */
      if(gVerbose)
      {
         printf("\nFinal Zones:\n");
         WriteTextOutput(zones, &numProts);
      }
      
      if(gDoOutput)
      {
         /*WRITE OUTPUT FILES FOR EACH PDB IN THE CORA ALIGNMENT*/
      }
   }
   else
   {
      Usage();
   }
   
   return(0);
}


/************************************************************************/
/*>BOOL ParseCmdLine(int argc, char **argv, char *ssapfile, 
                     char *pdbfile1, char *pdbfile2, char *outfile, 
                     char *outpdb1, char *outpdb2, REAL *dcut)
   ----------------------------------------------------------------
   Input:   int    argc         Argument count
            char   **argv       Argument array
   Output:  char   *ssapfile    Input SSAP file               
            char   *pdbfile1    First input PDB file
            char   *pdbfile2    Second input PDB file
            char   *outfile     Output listing file (or blank string)
            char   *outpdb1     Output first PDB file (or blank string)
            char   *outpdb2     Output second PDB file (or blank string)
            REAL   *dcut        Cutoff for defining core
   Returns: BOOL                Success?

   Parse the command line
   
   14.11.96 Original    By: ACRM
   06.12.96 Added -i
   23.01.97 Added -n
*/
BOOL ParseCmdLine(int argc, char **argv, char *corafile, REAL *dcut)
{
   argc--;
   argv++;
   corafile[0] = '\0';
   
   if(argc==0)
      return(FALSE);
   while(argc)
   {
      if (argv[0][0] == '-')
      {
         switch(argv[0][1])
         {
         case 'd':
            argc--;
            argv++;
            sscanf(argv[0],"%lf",dcut);
            break;
         case 'p':
            gDoOutput = TRUE;
            break;
         case 'v':
            gVerbose = TRUE;
            break;
         case 'i':
            gInitialCut = TRUE;
            break;
         case 'n':
            gDoRandomCoil = TRUE;
            break;
         default:
            return(FALSE);
            break;
         }
      }
      else
      {
         /* Check that there is  only 1 arguments left                  */
         if(argc < 1 || argc > 1)
            return(FALSE);
         
         /* Copy it                                                     */
         strcpy(corafile, argv[0]);
         
         return(TRUE);
      }
      argc--;
      argv++;
   }

   return(TRUE);
}

/************************************************************************/
/*>ZONE *ReadCORA(FILE *fp)
   ------------------------
   Read a CORA alignment file into a set of zones showing residue
   equivalences

   14.11.96 Original   By: ACRM 
   23.01.97 Added gDoRandomCoil checking; swapped the logic round for
            checking secondary structure matches to make this easier.
*/
Malign *ReadCORA(FILE *fp)
{
   /* Define a pointer to point to the first address in the Structure   */
   Malign *maln_ptr;
   /* Define a pointer to point to the first address in Malndata
      (loop through)
   */
   Malndata *d_ptr;          
   Protdata *p_ptr;  
   /* ZONE *zones; */
   
   int  count, count2, count3;
   char charbuf, buffer[1001], pdb[6];
   /*  char *p; */
   
   maln_ptr = new_Malign();
   while(fgets( buffer, MAXBUFF, fp ), buffer[0] == '#');
   
   /* read each line - maximum characters 100 - into the buffer string  */
   sscanf(buffer, "%d", &maln_ptr->procnt);
   
   /* copy the number of proteins to parse out of the function          */
   sprintf(maln_ptr->title, "%s", "alnfile");
   
   /* loop throught the number of aligned proteins and store the names  */
   for(count=0; count < maln_ptr->procnt; count++) 
   { 
      fscanf(fp, "%s ", maln_ptr->proname[count]);
   }
   /* store the number of lines in the file                             */
   fscanf(fp, "%d", &maln_ptr->length); 
   
   count = count2 = count3 = 0;
   /* Reserve space for the structure Malndata - it returns a pointer 
      which is made equal to the pointer malndata_ptr within the 
      structure Malign
   */
   (maln_ptr->malndata_ptr) = (Malndata*)malloc(sizeof(Malndata));
   TEST(maln_ptr->malndata_ptr);
      
   /* read through file saving data for selected chain                  */
   /* loop through every line in the file                               */
   while(count < maln_ptr->length) 
   { 	  
      /* reallocate space of size Malndata make it equal to the start
         of the structure
      */
      (maln_ptr->malndata_ptr) =
         (Malndata*)realloc(maln_ptr->malndata_ptr,
                            (count+1)*sizeof(Malndata)); 
      TEST(maln_ptr->malndata_ptr);
      
      /* move the pointer along to the next space by setting d_ptr      */
      d_ptr = maln_ptr->malndata_ptr + count; 
      
      /* store values                                                   */
      fscanf(fp, "%d %d %d ",
             &d_ptr->alnpos, &d_ptr->conpos, &d_ptr->proaln );
      
      /* Reserve space for the structure Protdata- it returns a pointer 
         which is made equal to the pointer malndata_ptr within the 
         structure d
      */
      (d_ptr->protdata_ptr)=(Protdata*)malloc(sizeof(Protdata));
      TEST(d_ptr->protdata_ptr);
      
      for(count2=0; count2 < maln_ptr->procnt; count2++) 
      {
         int  number;
         char insert;
         
         (d_ptr->protdata_ptr)=
            (Protdata*)realloc(d_ptr->protdata_ptr,
                               (count2+1)*sizeof(Protdata));
         TEST(d_ptr->protdata_ptr);
         
         p_ptr = d_ptr->protdata_ptr + count2;
         fscanf(fp, "%s %c %c ", pdb, &p_ptr->acid, &p_ptr->secstruct);
         insert = ' ';
         sscanf(pdb, "%d%c", &number, &insert);
         sprintf(p_ptr->pdb, "%d%c", number, insert);
      }
      
      fscanf(fp, "%c ", &d_ptr->consecstruc);
      while((charbuf = fgetc( fp )), charbuf != '\n' );
      count++;
   }
   
   return(maln_ptr);
}

/*************************************************************************
   FINDS THE ZONES IN CORA
*************************************************************************/
ZONE *calcZone(Malign *maln_ptr)
{
   int count, count2;
   Malndata *d_ptr;
   Malndata *d_temp_ptr;
   Malndata *d_temp_sec_ptr;
   Protdata *p_ptr;
   Protdata *p_temp_ptr;
   Protdata *p_temp_sec_ptr;
   
   ZONE *zones = NULL, *z;
   /*  int insertion = 0; */
   int strands = 0;
   int helix = 0; 
   int no_ss = 0;
   int startZone = 0;
   int endZone = 0;
   int protnum = 0;
   /*  int i = 0; */
   int first = 0;
   int second = 0;
   
   /* Outer loop through each alignment position                        */
   for(count2=0; count2 < maln_ptr->length; count2++) 
   {
      /* inner cycles through each protein                              */
      for( count=0; count < maln_ptr->procnt; count++)
      { 	
         d_ptr = maln_ptr->malndata_ptr + count2;	    
         p_ptr = d_ptr->protdata_ptr + count;
         
         if(p_ptr->secstruct == 'E')
            strands = 1;
         
         if(p_ptr->secstruct == 'H')
            helix = 1;

         if(p_ptr->secstruct == '0')
            no_ss = 1;
      }
      
      /* IF THERE ARE NO HELICES, COILS                                 */
      if ((strands == 1) && (no_ss == 0) && (helix == 0))
      {
         /* THIS IS THE START POSITION                                  */
         startZone = count2;
         /* while there are still only strands aligned and the end of
            the file is not reached
         */
         while((strands == 1) && (no_ss == 0) && (helix == 0) 
               && (count2 < maln_ptr->length))
         {
            /* go onto next alignment position                          */
            count2++;
            
            /* loop through all prots at that alignment position        */
            for (count=0; count < maln_ptr->procnt; count++)
            {
               d_ptr = maln_ptr->malndata_ptr + count2;	    
               p_ptr = d_ptr->protdata_ptr + count;
               if(p_ptr->secstruct == 'E')
                  strands = 1;

               if(p_ptr->secstruct == 'H')
                  helix = 1;

               if(p_ptr->secstruct == '0')
                  no_ss = 1;
            }
         }

         /* when this is not true any more mark the end of the zone
            minus one as the program increments one beyond the 
            zone to test it
         */
         endZone = count2;
         endZone = endZone-1;
      }
      /* IF THERE ARE NO STRANDS, COILS                                 */
      else if((strands == 0) && (no_ss == 0) && (helix == 1))
      {
         /* THIS IS THE START POSITION                                  */
         startZone = count2;
	 
         /* while there are still only strands aligned and the end of
            the file is not reached
         */
         while((strands == 0) && (no_ss == 0) && (helix == 1) 
               && (count2 < maln_ptr->length))
         {
            /* increment alignment position                             */
            count2++;
            
            /* loop through all prots at that alignment position        */
            for(count=0; count < maln_ptr->procnt; count++)
            {
               d_ptr = maln_ptr->malndata_ptr + count2;	    
               p_ptr = d_ptr->protdata_ptr + count;
               if(p_ptr->secstruct == 'E')
                  strands = 1;

               if(p_ptr->secstruct == 'H')
                  helix = 1;

               if(p_ptr->secstruct == '0')
                  no_ss = 1;
            }
         }

         /* when this is not true any more mark the end of the zone
            minus one as the program increments one beyond the 
            zone to test it
         */
         endZone = count2;
         endZone = endZone -1;
      }
      
      /* record the start and end residues into the zone array          */
      if(startZone)
      {
         if(zones==NULL)
         {
            INITPREV(zones,ZONE);
            z=zones;
         }
         else
         {
            ALLOCNEXTPREV(z,ZONE);
         }

         if(z==NULL)
         {
            FREELIST(zones,ZONE);
            return(NULL);
         }
         
         {
            int my_i;
            for(my_i=0; my_i<MAXMALNPNO; my_i++)
            {
               z->start[my_i] = (-1);
               z->end[my_i]   = (-1);
            }
         }
	 
         /* loop through the prots and record each resnum for start and
            end
         */
         for(protnum=0; protnum < maln_ptr->procnt; protnum++)
         {
            /* pointer to the start position                            */
            d_temp_ptr = maln_ptr->malndata_ptr + startZone;
            p_temp_ptr = d_temp_ptr->protdata_ptr + protnum;
            
            d_temp_sec_ptr = maln_ptr->malndata_ptr + endZone; 
            p_temp_sec_ptr = d_temp_sec_ptr->protdata_ptr + protnum;
            
            first = atoi(p_temp_ptr->pdb);
            second = atoi(p_temp_sec_ptr->pdb);
            
            z->start[protnum] = first;
            z->end[protnum] = second;
         } 
      }
      
      strands   = 0;
      helix     = 0;
      no_ss     = 0;
      startZone = 0;
      endZone   = 0;
   }
   
   return(zones);
}


/************************************************************************/
/*>BOOL DefineCore(FILE *outfp, PDB *pdb1, PDB *pdb2, ZONE *zones, 
                   REAL dcut)
  -----------------------------------------------------------------------
  Main routine to do core definition
  
  14.11.96 Original   By: ACRM
  06.12.96 Added handling of gInitialCut
*/
BOOL DefineCore(PDB **pdb, ZONE *zones, Malign *maln_ptr, REAL dcut)
{
   int  count    = 0,
        last     = 0,
        iter     = 0,
        protNum  = 0,
        numProts = 0,
        decrease = 0,
        natoms[MAXMALNPNO];
   REAL rm[3][3];
   PDB  *pdbca[MAXMALNPNO],
        **idx[MAXMALNPNO];
   
   /* Duplicate the PDB linked lists                                    */
   for(protNum = 0; protNum< maln_ptr->procnt; protNum++)
   {
      if((pdbca[protNum] = blDupePDB(pdb[protNum])) == NULL)
      {
         return(FALSE);
         for(decrease = protNum; decrease>1; decrease--)
         {
            FREELIST(pdbca[decrease],PDB);
         }
      }
   }
   
   /* Reduce to CA only                                                 */
   for(protNum = 0; protNum< maln_ptr->procnt; protNum++)
   {                    
      pdbca[protNum]= blSelectCaPDB(pdbca[protNum]);
      SetBValByZone(pdbca[protNum],zones,protNum);
   }

   count = CountCore(pdbca[0]);
   
   for(protNum = 0; protNum< maln_ptr->procnt; protNum++)
   {
      if((idx[protNum] = blIndexPDB(pdbca[protNum],
                                    &natoms[protNum]))==NULL)
      {
         for(decrease = protNum; decrease>1; decrease--)
         {
            FREELIST(pdbca[decrease],PDB);
            return(FALSE);
         }
      }
   }
   numProts = maln_ptr->procnt;
   
   if(gInitialCut)
   {
      for(protNum=1; protNum<numProts; protNum++)
      {
         FitCaPDBBFlag(pdbca[0], pdbca[protNum], rm);
      }
      
      if(!DoCut(idx, natoms, numProts, zones, dcut*dcut))
         return(FALSE);
      
      count = CountCore(pdbca[0]);
      
      if(gVerbose)
      {
         printf("\nCore after removing residues > 3.0A:\n");
         WriteTextOutput(zones, &numProts);
      }
   }
   
   iter=0;
   
   while(last != count)
   {
      for(protNum=1; protNum<numProts; protNum++)
      {
         FitCaPDBBFlag(pdbca[0], pdbca[protNum], rm);
      }
      
      UpdateBValues(idx, natoms, numProts, zones, dcut*dcut);
      last = count;
      count = CountCore(pdbca[0]);
      if(++iter > MAXITER)
      {
         fprintf(stderr,"Warning: Maximum number of iterations \
(%d) exceeded!\n",MAXITER);
         break;
      }
   }
   
   for(protNum = 0; protNum < maln_ptr->procnt; protNum++)
   {
      free(idx[protNum]);
   }
   
   return(TRUE);
}


/************************************************************************/
/*>BOOL DoCut(PDB **idx1, int natom1, PDB **idx2, int natom2,
  ZONE *zones, REAL cutsq)
  ----------------------------------------------------------
  Performs the initial cut of pairs which deviate by >3.0A
  
  06.12.96 Original   By: ACRM
*/
BOOL DoCut(PDB **idx[MAXMALNPNO], int *natoms, int nstruc,
           ZONE *zones, REAL cutsq)
{
   ZONE *z, *zend, *znext;
   int  i,
        snum, snum2, snum3,
        starts[MAXMALNPNO], ends[MAXMALNPNO], offsets[MAXMALNPNO];
   BOOL split,
        ok,
        lastIter;
   
   for(z=zones; z!=NULL; NEXT(z))
   {
      /* Look for the start of the zone                                 */
      for(snum=0; snum<nstruc; snum++)
      {
         for(starts[snum]=0; starts[snum]<natoms[snum]; starts[snum]++)
         {
            if(idx[snum][starts[snum]]->resnum == z->start[snum])
               break;
         }
      }
      
      /* Look for the end of the zone                                   */
      for(snum=0; snum<nstruc; snum++)
      {
         for(ends[snum]=0; ends[snum]<natoms[snum]; ends[snum]++)
         {
            if(idx[snum][ends[snum]]->resnum == z->end[snum])
               break;
         }
      }
      
      /* Step through the zone to see if we are within the cutoff       */
      split = FALSE;
      for(snum=0; snum<nstruc; snum++)
      {
         offsets[snum] = starts[snum];
      }
      
      lastIter = FALSE;
      for(;;)
      {
         /* Check the offsets are all in range                          */
         for(snum=0; snum<nstruc; snum++)
         {
            if(offsets[snum] > ends[snum])
            {
               lastIter = TRUE;
               break;
            }
         }
         if(lastIter)
            break;
         
         /* Check all the pairwise distances are in range               */
         for(snum=0; snum<nstruc-1; snum++)
         {
            for(snum2=snum+1; snum2<nstruc; snum2++)
            {
               /* If a pair is out of range, then set the split flag and
                  set the bvalues of all the members to zero
               */
               if(DISTSQ(idx[snum][offsets[snum]],
                         idx[snum2][offsets[snum2]]) > cutsq)
               {
                  split = TRUE;
                  for(snum3=0; snum3<nstruc; snum3++)
                  {
                     idx[snum3][offsets[snum3]]->bval = (REAL)0.0;
                  }
               }
            }
         }
         
         /* Increment the offsets                                       */
         for(snum=0; snum<nstruc; snum++)
         {
            (offsets[snum])++;
         }
      }
      
      /* See if there were any bits out of range which need the zones
         to be split
      */
      if(split)
      {
         /* See if the whole zone was out of range                      */
         ok = FALSE;
         lastIter = FALSE;
         for(snum=0; snum<nstruc; snum++)
         {
            for(i=starts[snum]; i<=ends[snum]; i++)
            {
               if(idx[snum][i]->bval == (REAL)10.0)
               {
                  ok = TRUE;
                  lastIter = TRUE;
                  break;
               }
            }
            if(lastIter)
               break;
         }
         
         /* If the zone contains no pairs within 3.0A, mark it for
            deletion by setting all values to -9999
         */
         if(!ok)
         {
            for(snum=0; snum<nstruc; snum++)
            {
               z->start[snum] = z->end[snum] = (-9999);
            }
         }
         else
         {
            /* There were some parts which are still required, but it's
               been modified, so we need to split the zones up
            */
            
            /* First see if we've lost residues from the start of the
               zone
            */
            for(;;)
            {
               /* Start off assuming that all the pairs are OK at the 
                  start of the zone and we are not removing anything.
                  See if any of the pairs were to be removed and if so
                  this isn't the last iteration
               */
               lastIter = TRUE;
               for(snum=0; snum<nstruc; snum++)
               {
                  if(idx[snum][starts[snum]]->bval == (REAL)0.0)
                  {
                     lastIter = FALSE;
                     break;
                  }
               }
               
               /* If it's the last iteration, so nothing more to be
                  removed from the start of the zone, then break out
               */
               if(lastIter)
               {
                  break;
               }
               else
               {
                  /* We are removing this residue group from the zone
                     so bump the start pointers and update the zone
                  */
                  for(snum=0; snum<nstruc; snum++)
                  {
                     (starts[snum])++;
                     z->start[snum] = idx[snum][starts[snum]]->resnum;
                  }
               }
            }
            
            /* Now remove residues from the end of the zone in the same
               way
            */
            for(;;)
            {
               /* Start off assuming that all the pairs are OK at the 
                  start of the zone and we are not removing anything.
                  See if any of the pairs were to be removed and if so
                  this isn't the last iteration
               */
               lastIter = TRUE;
               for(snum=0; snum<nstruc; snum++)
               {
                  if(idx[snum][ends[snum]]->bval == (REAL)0.0)
                  {
                     lastIter = FALSE;
                     break;
                  }
               }
               
               /* If it's the last iteration, so nothing more to be
                  removed from the start of the zone, then break out
               */
               if(lastIter)
               {
                  break;
               }
               else
               {
                  /* We are removing this residue group from the zone
                     so bump the start pointers and update the zone
                  */
                  for(snum=0; snum<nstruc; snum++)
                  {
                     (ends[snum])--;
                     z->end[snum] = idx[snum][ends[snum]]->resnum;
                  }
               }
            }
            
            /* See if the new zone is split                             */
            split = FALSE;
            for(snum=0; snum<nstruc; snum++)
            {
               for(i=starts[snum]; i<=ends[snum]; i++)
               {
                  if(idx[snum][i]->bval == (REAL)0.0)
                  {
                     split = TRUE;
                     break;
                  }
               }
               if(split)
                  break;
            }
            
            while(split)
            {
               /* Move the last sub-zone into a separate zone           */
               znext = z->next;
               INIT(zend, ZONE);
               if(zend==NULL)
               {
                  fprintf(stderr,"No memory for new zones!\n");
                  return(FALSE);
               }
               z->next     = zend;
               zend->prev  = z;
               zend->next  = znext;
               znext->prev = zend;
               
               for(snum=0; snum<nstruc; snum++)
               {
                  zend->end[snum] = z->end[snum];
               }
               
               /* Step back through the zone to find the start of this
                  subzone
               */
               lastIter = FALSE;
               for(;;)
               {
                  for(snum=0; snum<nstruc; snum++)
                  {
                     if(idx[snum][ends[snum]]->bval != (REAL)10.0)
                     {
                        lastIter = TRUE;
                        break;
                     }
                  }
                  if(lastIter)
                     break;
                  
                  for(snum=0; snum<nstruc; snum++)
                  {
                     (ends[snum])--;
                     zend->start[snum] = idx[snum][ends[snum]]->resnum;
                  }
               }
               
               /* Now step back to the end of the previous subzone      */
               lastIter = FALSE;
               for(;;)
               {
                  for(snum=0; snum<nstruc; snum++)
                  {
                     if(idx[snum][ends[snum]]->bval != (REAL)0.0)
                     {
                        lastIter = TRUE;
                        break;
                     }
                  }
                  if(lastIter)
                     break;
                  
                  for(snum=0; snum<nstruc; snum++)
                  {
                     (ends[snum])--;
                     z->end[snum] = idx[snum][ends[snum]]->resnum;
                  }
               }
               
               /* Test again to see if it's split                       */
               split = FALSE;
               
               for(snum=0; snum<nstruc; snum++)
               {
                  for(i=starts[snum]; i<=ends[snum]; i++)
                  {
                     if(idx[snum][i]->bval == (REAL)0.0)
                     {
                        split = TRUE;
                        break;
                     }
                  }
                  if(split)
                     break;
               }
            }  /* While the zone is split                               */
         }  /* If this zone was OK                                      */
      }  /* If the zone contained any splits                            */
   }  /* For each zone                                                  */
   
   return(TRUE);
}

/************************************************************************/
/*>void UpdateBValues(PDB **idx[MAXSTRUC], int *natoms, int nstruc,
                      ZONE *zones, REAL cutsq)
   ------------------------------------------------------------------
   Update the B-values and the current zones by extending out from
   the secondary structure regions
  
   14.11.96 Original   By: ACRM
   13.03.97 Now checks that residues are not already in a zone before
            adding them to the current zone. Fixes a problem at the
            zone-merge stage where the merged zones could end up with
            different numbers of residues.
   08.05.02 Generalized to work with multiple structures
*/
void UpdateBValues(PDB **idx[MAXMALNPNO], int *natoms, int nstruc,
                   ZONE *zones, REAL cutsq)
{
   ZONE *z;
   int  snum, offset[MAXMALNPNO];
   BOOL lastIter;
   
   for(z=zones; z!=NULL; NEXT(z))
   {
      if(z->start[0] < -9998)
         continue;
      
      /* Look for the start of the zone - store this in offset[] for
         each protein
      */
      for(snum=0; snum<nstruc; snum++)
      {
         for(offset[snum]=0; offset[snum]<natoms[snum]; offset[snum]++)
         {
            if(idx[snum][offset[snum]]->resnum == z->start[snum])
               break;
         }
      }
      
      /* Step back from the start of a zone, checking every pair of
         atoms is within the cutoff. Keep going back until an atom pair
         is too far apart or we reach the beginning of one of the proteins
      */
      
      /* Loop forever                                                   */
      lastIter = FALSE;
      for(;;)
      {
         /* Move the protein pointers back one residue, checking that we 
            haven't come off the start of any of the structures.
         */
         for(snum=0; snum<nstruc; snum++)
         {
            if(--offset[snum] < 0)
            {
               lastIter = TRUE;
               /* Note that we don't break out here as we need all the
                  pointers to be in sync
               */
            }
            
         }
         if(lastIter)
            break;
         
         /* End the loop if the BVal for any of the structures is >5.0  */
         for(snum=0; snum<nstruc; snum++)
         {
            if(idx[snum][offset[snum]]->bval > (REAL)5.0)
            {
               lastIter = TRUE;
               break;
            }
         }
         if(lastIter)
            break;
         
         /* Now look at each pair of structures, if the distance between
            the atoms is > cutsq then end the loop
         */
         for(snum=0; snum<nstruc-1; snum++)
         {
            int snum2;
            for(snum2=snum+1; snum2<nstruc; snum2++)
            {
               if((DISTSQ(idx[snum][offset[snum]],
                          idx[snum2][offset[snum2]])) > cutsq)
               {
                  lastIter = TRUE;
                  snum = nstruc; /* Force break from outer loop         */
                  break;
               }
            }
         }
         if(lastIter)
            break;
         
         /* Everything was OK, so mark these residues as part of a zone */
         for(snum=0; snum<nstruc; snum++)
         {
            idx[snum][offset[snum]]->bval = (REAL)10.0;
         }
      }  /* End of loop back from start of zone                         */
      
      /* Bump the offsets back up by one so we are pointing to the first
         residue that was within the OK zone and update the zone starts
      */
      for(snum=0; snum<nstruc; snum++)
      {
         (offset[snum])++;
         z->start[snum] = idx[snum][offset[snum]]->resnum;
      }
      
      /* ---- Repeat the above procedure but for the ends of zones ---- */
      
      /* Look for the end of the zone - store this in offset[] for each
         protein
      */
      for(snum=0; snum<nstruc; snum++)
      {
         for(offset[snum]=0; offset[snum]<natoms[snum]; offset[snum]++)
         {
            if(idx[snum][offset[snum]]->resnum == z->end[snum])
               break;
         }
      }
      
      /* Step forward from the end of zone, checking every pair of
         atoms is within the cutoff. Keep going forward until an atom 
         pair is too far apart or we reach the end of one of the proteins
      */
      
      /* Loop forever                                                   */
      lastIter = FALSE;
      for(;;)
      {
         /* Move the protein pointers forward one residue, checking that
            we haven't come off the end of any of the structures.
         */
         for(snum=0; snum<nstruc; snum++)
         {
            if(++offset[snum] >= natoms[snum])
            {
               lastIter = TRUE;
               /* Note that we don't break out here as we need all the
                  pointers to be in sync
               */
            }
            
         }
         if(lastIter)
            break;
         
         /* End the loop if the BVal for any of the structures is >5.0  */
         for(snum=0; snum<nstruc; snum++)
         {
            if(idx[snum][offset[snum]]->bval > (REAL)5.0)
            {
               lastIter = TRUE;
               break;
            }
         }
         if(lastIter)
            break;
         
         /* Now look at each pair of structures, if the distance between
            the atoms is > cutsq then end the loop
         */
         for(snum=0; snum<nstruc-1; snum++)
         {
            int snum2;
            for(snum2=snum+1; snum2<nstruc; snum2++)
            {
               if((DISTSQ(idx[snum][offset[snum]],
                          idx[snum2][offset[snum2]])) > cutsq)
               {
                  lastIter = TRUE;
                  snum = nstruc; /* Force break from outer loop         */
                  break;
               }
            }
         }
         if(lastIter)
            break;
         
         /* Everything was OK, so mark these residues as part of a zone */
         for(snum=0; snum<nstruc; snum++)
         {
            idx[snum][offset[snum]]->bval = (REAL)10.0;
         }
      }  /* End of loop forwards from end of zone                       */
      
      /* Bump the offsets back down by one so we are pointing to the last
         residue that was within the OK zone and update the zone ends
      */
      for(snum=0; snum<nstruc; snum++)
      {
         (offset[snum])--;
         z->end[snum] = idx[snum][offset[snum]]->resnum;
      }
   }  /* Continue with the next zone                                    */
}


/************************************************************************/
/*>void SetBValByZone(PDB *pdb, ZONE *zones, int which)
  ----------------------------------------------------
  Set B-values to 10 if in the zones, otherwise to 0.0
  
  14.11.96 Original   By: ACRM
*/
void SetBValByZone(PDB *pdb, ZONE *zones, int protNum)
{
   PDB *p;
   ZONE *z;
   
   for(p=pdb; p!=NULL; NEXT(p))
      p->bval = (REAL)0.0;
   
   for(p=pdb; p!=NULL; NEXT(p))
   {
      for(z=zones; z!=NULL; NEXT(z))
      {
         if((p->resnum >= z->start[protNum]) &&
            (p->resnum <= z->end[protNum]))
            p->bval = (REAL)10.0;
      }
   }
}


/************************************************************************/
/*>BOOL FitCaPDBBFlag(PDB *ref_pdb, PDB *fit_pdb, REAL rm[3][3])
  -------------------------------------------------------------
  Input:   PDB  *ref_pdb     Reference PDB linked list
  I/O:     PDB  *fit_pdb     Mobile PDB linked list
  Output:  REAL rm[3][3]     Rotation matrix (May be input as NULL).
  Returns: BOOL              Success
  
  Fits two PDB linked lists using only the CA atoms of residues
  flagged in the BVal column with a non-zero BValue (atoms not to be 
  fitted have BVal set to 0).
  
  Actually fits fit_pdb onto ref_pdb and also returns the rotation 
  matrix. This may be NULL if these data are not required.
  
  14.11.96 Original based on FitCaPDB()   By: ACRM
*/
BOOL FitCaPDBBFlag(PDB *ref_pdb, PDB *fit_pdb, REAL rm[3][3])
{
   REAL  RotMat[3][3];
   COOR  *ref_coor   = NULL,
         *fit_coor   = NULL;
   VEC3F ref_ca_CofG,
         fit_ca_CofG,
         tvect;
   int   NCoor       = 0,
         i, j;
   BOOL  RetVal;
   PDB   *ref_ca_pdb = NULL,
         *fit_ca_pdb = NULL;
   
   /* First extract only the CA atoms of residues where BVal > 0.5      */
   if((ref_ca_pdb = DupeCAByBVal(ref_pdb))==NULL)
      RetVal = FALSE;

   if((fit_ca_pdb = DupeCAByBVal(fit_pdb))==NULL)
      RetVal = FALSE;
   
   /* If we succeeded in building our CA PDB linked lists...            */
   if(RetVal)
   {
      /* Get the CofG of the CA structures and the original mobile      */
      blGetCofGPDB(ref_ca_pdb, &ref_ca_CofG);
      blGetCofGPDB(fit_ca_pdb, &fit_ca_CofG);
      
      /* Move them both to the origin                                   */
      blOriginPDB(ref_ca_pdb);
      blOriginPDB(fit_ca_pdb);
      /*
      for(p=fit_ca_pdb;  p!=NULL; NEXT(p))
      {
         printf("%f  \n", p->x);
      }
      */

      /* Create coordinate arrays, checking numbers match               */
      NCoor = blGetPDBCoor(ref_ca_pdb, &ref_coor);
      if(blGetPDBCoor(fit_ca_pdb, &fit_coor) != NCoor)
      {
         RetVal = FALSE;
      }
      else
      {
         /* Can't fit with fewer than 3 coordinates                     */
         if(NCoor < 3)
         {
            RetVal = FALSE;
         }
         else
         {
            /* Everything OK, go ahead with the fitting                 */
            if(!blMatfit(ref_coor,fit_coor,RotMat,NCoor,NULL,FALSE))
            {
               RetVal = FALSE;
            }
            else
            {
               /* Apply the operations to the true coordinates          */
               tvect.x = (-fit_ca_CofG.x);
               tvect.y = (-fit_ca_CofG.y);
               tvect.z = (-fit_ca_CofG.z);
               blTranslatePDB(fit_pdb, tvect);
               blApplyMatrixPDB(fit_pdb, RotMat);
               blTranslatePDB(fit_pdb, ref_ca_CofG);
            }
         }
      }
   }
   
   /* Free the coordinate arrays and CA PDB linked lists                */
   if(ref_coor)   free(ref_coor);
   if(fit_coor)   free(fit_coor);
   if(ref_ca_pdb) FREELIST(ref_ca_pdb, PDB);
   if(fit_ca_pdb) FREELIST(fit_ca_pdb, PDB);
   
   /* Fill in the rotation matrix for output, if required               */
   if(RetVal && (rm!=NULL))
   {
      for(i=0; i<3; i++)
         for(j=0; j<3; j++)
            rm[i][j] = RotMat[i][j];
   }
   
   return(RetVal);
}


/************************************************************************/
/*>int CountCore(PDB *pdb)
  -----------------------
  Count how many residues are in the core regions
  
  14.11.96 Original   By: ACRM
*/
int CountCore(PDB *pdb)
{
   PDB *p;
   int count = 0;
   
   for(p=pdb; p!=NULL; NEXT(p))
      if(p->bval > (REAL)0.0)
         count++;
   
   return(count);
}


/************************************************************************/
/*>PDB *DupeCAByBVal(PDB *pdb)
  ---------------------------
  Extract the CAs with BVal > 0.0 from a PDB linked list into a new
  linked list
  
  14.11.96 Original   By: ACRM
*/
PDB *DupeCAByBVal(PDB *pdb)
{
   PDB *out = NULL, 
       *p, *q;
   
   for(p=pdb; p!=NULL; NEXT(p))
   {
      if(!strncmp(p->atnam,"CA  ",3) && (p->bval > (REAL)0.0))
      {
         if(out==NULL)
         {
            INIT(out,PDB);
            q=out;
         }
         else
         {
            ALLOCNEXT(q,PDB);
         }
         if(q==NULL)
         {
            FREELIST(out,PDB);
            return(NULL);
         }
         
         blCopyPDB(q,p);
      }
   }
   
   return(out);
}

Malign *new_Malign(void)
{
   Malign *maln;
   ALLOC( maln, Malign );
   TEST(maln);
   clear_Malign(maln);
   return(maln);
}

void clear_Malign(Malign *m)
{
   int i;
   
   m->procnt = 0;
   m->length = 0;
   *(m->title) = '\0';
   
   for(i=0; i<MAXMALNPNO; i++) 
   {
      *(m->proname[i]) = '\0';
   }
}


/************************************************************************/
/*>void WriteTextOutput(FILE *fp, ZONE *zones)
  -------------------------------------------
  Writes the zones out in text format
  
  14.11.96 Original   By: ACRM
  06.12.96 Added check that zones have not been blanked out
*/
void WriteTextOutput(ZONE *zones, int *numProts)
{
   ZONE *z;
   int  i = 0;
   
   for(z=zones; z!=NULL; NEXT(z))
   {
      if(z->start[0] > -9999)
      {
	 for(i = 0; i< *numProts; i++)
         {
            printf("%4d -%4d ",z->start[i],z->end[i]);
            if (i<*numProts)
            {
               printf(": ");
            }
         }
	 printf("\n");
      }
   }
}


/************************************************************************/
/*>void Usage(void)
  ----------------
  Prints a usage message
  
  14.11.96 Original   By: ACRM
  06.12.96 V1.1
  23.01.97 V1.2
  05.11.25 V1.8
*/
void Usage(void)
{
   fprintf(stderr,"\nFindCore V1.8 (c) 1996-2025, Prof. Andrew C.R. \
Martin, UCL.\n");
   fprintf(stderr,"Modifications for Cora by Gabby Marsden (nee Reeves) \
           1999-2002\n");
   
   fprintf(stderr,"\nUsage: findcore [-p out1.pdb] [-q out2.pdb] [-d \
dcut] [-v] [-i]\n");
   fprintf(stderr,"                ssapfile in1.pdb in2.pdb \
[output.lis]\n");
   fprintf(stderr,"       -p       Write in1.pdb with core flagged in \
B-value column\n");
   fprintf(stderr,"       -q       Write in2.pdb with core flagged in \
B-value column\n");
   fprintf(stderr,"       -d       Specify distance cutoff for defining \
core [%f]\n",DEFAULT_CUT);
   fprintf(stderr,"       -v       Verbose mode; shows intermediate \
zones\n");
   fprintf(stderr,"       -i       Initial cut. Do an initial fit and \
remove any pairs\n");
   fprintf(stderr,"                with > spcified cutoff\n");
   fprintf(stderr,"       -n       Include non-E/H regions which match \
in the initial\n");
   fprintf(stderr,"                definition of core zones\n");
   fprintf(stderr,"       ssapfile A vertical alignment file from \
SSAP\n");
   
   fprintf(stderr,"\nFindCora defines a protein structurally conserved \
core according to the\n");
   fprintf(stderr,"method of Chothia (1996). Structurally equivalent \
secondary structure\n");
   fprintf(stderr,"regions (defined by SSAP) are fitted. Residues are \
added at the end\n");
   fprintf(stderr,"of each region while the CA-deviation is less than \
dcut and the fitting\n");
   fprintf(stderr,"is repeated. This iterates until no additional \
residues are added.\n\n");
   fprintf(stderr,"Note that the program will not handle chain names or \
insertion codes in\n");
   fprintf(stderr,"the PDB files. This is a limitation imposed by \
SSAP.\n\n");
   
   fprintf(stderr,"The PDB files should be given in the same order as \
the columns appear\n");
   fprintf(stderr,"in the SSAP file.\n\n");
}


/************************************************************************/
/*>ZONE *MergeZones(ZONE *zones, int numProts)
  -------------------------------------------
  Merges zones in the zone linked list of they overlap and removes zones
  which are subsets of other zones
  
  14.11.96 Original   By: ACRM
  06.12.96 Removes zones marked for deletion
  13.03.97 Added merging of abutting zones (required since change to
           zone creation where residues not added to a zone if they
           are already in another zone stops them from being subsets)
  26.06.02 Generalized for multiple structures
           Fixed bug in deleting zones
*/
ZONE *MergeZones(ZONE *zones, int numProts)
{
   ZONE *z, *zt;
   BOOL finished = FALSE,
        doit;
   int  i;
   
   /* Remove null zones                                                 */
   while(zones->start[0] < -9998)
      zones = zones->next;
   for(z=zones;z!=NULL;NEXT(z))
   {
      if(z->start[0] < -9998)
      {
         ZONE *prev;
         
         z->prev->next = z->next;
         if(z->next != NULL)
            z->next->prev = z->prev;
         prev = z->prev;
         free(z);
         z = prev;
      }
   }
   
   /* First merge overlapping zones                                     */
   while(!finished)
   {
      finished = TRUE;
      for(z=zones; z->next!=NULL; NEXT(z))
      {
         /* For each structure, see if the end of a zone overlaps the 
            start of a following zone
         */
         for(i=0, doit=FALSE; i<numProts; i++)
         {
            if(z->end[i] >= z->next->start[i])
            {
               doit = TRUE;
               break;
            }
         }
         
         if(doit)  /* There was a zone overlap                          */
         {
            /* For each structure, check that the end of the zone doesn't
               already equal the end of the following zone. If there
               weren't any such cases, then we will merge zones.
            */
            for(i=0, doit=TRUE; i<numProts; i++)
            {
               if(z->end[i] == z->next->end[i])
               {
                  doit = FALSE;
                  break;
               }
            }
            
            if(doit) /* The zone didn't already span the following zone */
            {
               finished = FALSE;
               for(i=0; i<numProts; i++)
               {
                  z->end[i] = z->next->end[i];
               }
            }
         }
      }
   }
   
   /* Now remove redundant zones                                        */
   finished = FALSE;
   while(!finished)
   {
      finished = TRUE;
      for(z=zones; z!=NULL; NEXT(z))
      {
         if(SubsetZone(z, zones))
         {
            if(z==zones)
            {
               zones=zones->next;
               zones->prev=NULL;
            }
            else
            {
               z->prev->next = z->next;
               if(z->next != NULL)
                  z->next->prev = z->prev;
            }
            zt=z->prev;
            free(z);
            z=zt;
            
            finished = FALSE;
            break;
         }
      }
   }
   
   /* Now merge abutting zones                                          */
   finished = FALSE;
   while(!finished)
   {
      finished = TRUE;
      for(z=zones; z->next!=NULL; NEXT(z))
      {
         /* See if this zone abutts the next zone. Assume that it does,
            then check through each structure to see if there is a 
            mismatch
         */
         for(i=0, doit=TRUE; i<numProts; i++)
         {
            if(z->end[i]+1 != z->next->start[i])
            {
               doit = FALSE;
               break;
            }
         }
         
         if(doit) /* If the zones were abutting                         */
         {
            /* Reset start of second zone to start of first zone        */
            for(i=0; i<numProts; i++)
            {
               z->next->start[i] = z->start[i];
            }
            
            /* Remove z from the linked list                            */
            if(z==zones)
            {
               zones=zones->next;
               zones->prev=NULL;
            }
            else
            {
               if(z->prev != NULL) z->prev->next = z->next;
               if(z->next != NULL)
                  z->next->prev = z->prev;
            }
            zt=z->prev;
            free(z);
            z=zt;
            
            finished = FALSE;
         }
      }
   }
   
   return(zones);
}

/************************************************************************/
/*>BOOL SubsetZone(ZONE *z, ZONE *zones)
  -------------------------------------
  Sees whether a zones is a zubset of another zone
  
  14.11.96 Original   By: ACRM
*/
BOOL SubsetZone(ZONE *z, ZONE *zones)
{
   ZONE *z1;
   
   for(z1=zones; z1!=NULL; NEXT(z1))
   {
      if(z1 != z)
      {
         if((z->start[0] >= z1->start[0]) &&
            (z->end[0]   <= z1->end[0]))
            return(TRUE);
      }
   }
   
   return(FALSE);
}
