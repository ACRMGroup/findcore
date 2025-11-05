/*************************************************************************

   Program:    findcore
   File:       findcore.c
   
   Version:    V1.5
   Date:       05.11.25
   Function:   Find core from 2 structures given the SSAP alignment
               file as a staring point
   
   Copyright:  (c) Prof. Andrew C. R. Martin, UCL 1996-2025
   Author:     Prof. Andrew C. R. Martin
   Address:    Biomolecular Structure & Modelling Unit,
               Department of Biochemistry & Molecular Biology,
               University College,
               Gower Street,
               London.
               WC1E 6BT.
   EMail:      andrew@bioinf.org.uk
               
**************************************************************************

   This program is licensed under the GPL.

**************************************************************************

   Description:
   ============

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
   V1.4  26.06.02 Fixed bug in freeing zones in MergeZone()
   V1.5  05.11.25 Updated for new BiopLib

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
#define DEFAULT_CUT ((REAL)3.0)
typedef struct _zone
{
   struct _zone *next, *prev;
   int start[2], 
       end[2];
}  ZONE;

/************************************************************************/
/* Globals
*/
BOOL gVerbose      = FALSE,
     gInitialCut   = FALSE,
     gDoRandomCoil = FALSE;

/************************************************************************/
/* Prototypes
*/
int main(int argc, char **argv);
BOOL ParseCmdLine(int argc, char **argv, char *ssapfile, char *pdbfile1,
                  char *pdbfile2, char *outfile, char *outpdb1, 
                  char *outpdb2, REAL *dcut);
int strlen_nospace(char *str);
ZONE *ReadSSAP(FILE *fp);
BOOL DefineCore(FILE *outfp, PDB *pdb1, PDB *pdb2, ZONE *zones, REAL dcut);
void UpdateBValues(PDB **idx1, int natom1, PDB **idx2, int natom2,
                   ZONE *zones, REAL cutsq);
void SetBValByZone(PDB *pdb, ZONE *zones, int which);
BOOL FitCaPDBBFlag(PDB *ref_pdb, PDB *fit_pdb, REAL rm[3][3]);
int CountCore(PDB *pdb);
PDB *DupeCAByBVal(PDB *pdb);
void Usage(void);
void WriteTextOutput(FILE *fp, ZONE *zones);
BOOL SubsetZone(ZONE *z, ZONE *zones);
ZONE *MergeZones(ZONE *zones);
BOOL DoCut(PDB **idx1, int natom1, PDB **idx2, int natom2,
           ZONE *zones, REAL cutsq);



/************************************************************************/
/*>int main(int argc, char **argv)
   -------------------------------
   Main program for core defining

   14.11.96 Original   By: ACRM
*/
int main(int argc, char **argv)
{
   char ssapfile[MAXBUFF],
        pdbfile1[MAXBUFF],
        pdbfile2[MAXBUFF],
        outfile[MAXBUFF],
        outpdb1[MAXBUFF],
        outpdb2[MAXBUFF];
   FILE *ssapfp,
        *pdb1fp,
        *pdb2fp,
        *outfp = stdout;
   REAL dcut = DEFAULT_CUT;
   PDB  *pdb1,
        *pdb2;
   int  natoms;
   ZONE *zones;

   if(ParseCmdLine(argc, argv, ssapfile,pdbfile1,pdbfile2,outfile,
                   outpdb1,outpdb2,&dcut))
   {
      /* Open files                                                     */
      if((ssapfp=fopen(ssapfile,"r"))==NULL)
      {
         fprintf(stderr,"Unable to open %s for reading\n",ssapfile);
         return(1);
      }
      if((pdb1fp=fopen(pdbfile1,"r"))==NULL)
      {
         fprintf(stderr,"Unable to open %s for reading\n",pdbfile1);
         return(1);
      }
      if((pdb2fp=fopen(pdbfile2,"r"))==NULL)
      {
         fprintf(stderr,"Unable to open %s for reading\n",pdbfile2);
         return(1);
      }
      if(outfile[0])
      {
         if((outfp = fopen(outfile,"w"))==NULL)
         {
            fprintf(stderr,"Unable to open %s for writing\n",outfile);
            return(1);
         }
      }

      /* Read data from the files                                       */
      if((pdb1 = blReadPDB(pdb1fp,&natoms)) == NULL)
      {
         fprintf(stderr,"No atoms read from PDB file: %s\n",pdbfile1);
         return(1);
      }
      fclose(pdb1fp);
      if((pdb2 = blReadPDB(pdb2fp,&natoms)) == NULL)
      {
         fprintf(stderr,"No atoms read from PDB file: %s\n",pdbfile2);
         return(1);
      }
      fclose(pdb2fp);
      if((zones = ReadSSAP(ssapfp))==NULL)
      {
         fprintf(stderr,"No zones read from SSAP file: %s\n",ssapfile);
         return(1);
      }
      fclose(ssapfp);

      /* Print the current zones if required                            */
      if(gVerbose)
      {
         fprintf(outfp,"SSAP Zones:\n");
         WriteTextOutput(outfp, zones);
      }

      /* Now call the routine to do the core definition                 */
      DefineCore(outfp, pdb1, pdb2, zones, dcut);
      
      if(gVerbose)
      {
         fprintf(outfp,"\nCore before zone merging:\n");
         WriteTextOutput(outfp, zones);
      }
      
      /* Now remove any zones which are subsets of other zones and merge
         overlapping zones
      */
      zones = MergeZones(zones);
      
      /* Finally write the output file which lists residues in the 
         structural core and optionally write PDB files with the cores
         flagged
      */
      if(gVerbose)
         fprintf(outfp,"\nFinal Zones:\n");
      WriteTextOutput(outfp, zones);

      if(outpdb1[0])
      {
         if((pdb1fp=fopen(outpdb1,"w"))==NULL)
         {
            fprintf(stderr,"Unable to open %s for writing\n",outpdb1);
            return(1);
         }
         SetBValByZone(pdb1, zones, 0);
         blWritePDB(pdb1fp, pdb1);
      }
      if(outpdb2[0])
      {
         if((pdb2fp=fopen(outpdb2,"w"))==NULL)
         {
            fprintf(stderr,"Unable to open %s for writing\n",outpdb2);
            return(1);
         }
         SetBValByZone(pdb2, zones, 1);
         blWritePDB(pdb2fp, pdb2);
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
BOOL ParseCmdLine(int argc, char **argv, char *ssapfile, char *pdbfile1,
                  char *pdbfile2, char *outfile, char *outpdb1, 
                  char *outpdb2, REAL *dcut)
{
   argc--;
   argv++;

   ssapfile[0] = pdbfile1[0] = pdbfile2[0] = 
      outfile[0] = outpdb1[0] = outpdb2[0] = '\0';

   if(argc==0)
      return(FALSE);
   
   while(argc)
   {
      if(argv[0][0] == '-')
      {
         switch(argv[0][1])
         {
         case 'd':
            argc--;
            argv++;
            sscanf(argv[0],"%lf",dcut);
            break;
         case 'p':
            argc--;
            argv++;
            strcpy(outpdb1, argv[0]);
            break;
         case 'q':
            argc--;
            argv++;
            strcpy(outpdb2, argv[0]);
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
         /* Check that there are only 3 or 4 arguments left             */
         if(argc < 3 || argc > 4)
            return(FALSE);
         
         /* Copy the first three                                        */
         strcpy(ssapfile, argv[0]);
         strcpy(pdbfile1, argv[1]);
         strcpy(pdbfile2, argv[2]);
         
         /* If there's another, copy it to outfile                      */
         argc -= 3;
         argv += 3;
         if(argc)
            strcpy(outfile, argv[0]);
            
         return(TRUE);
      }
      argc--;
      argv++;
   }
   
   return(TRUE);
}


/************************************************************************/
/*>int strlen_nospace(char *str)
   -----------------------------
   Find the length of a string without any leading or trailing spaces

   14.11.96 Original   By: ACRM
*/
int strlen_nospace(char *str)
{
   char *s,
        *e;
   int  str_len;
   
   for(s=str; *s==' ' || *s=='\t'; s++);
   str_len = strlen(s);
   for(e=s+str_len-1; *e==' ' || *e=='\t'; e--)
      str_len--;
   return(str_len);
}

   
/************************************************************************/
/*>ZONE *ReadSSAP(FILE *fp)
   ------------------------
   Read a SSAP alignment file into a set of zones showing residue
   equivalences

   14.11.96 Original   By: ACRM
   23.01.97 Added gDoRandomCoil checking; swapped the logic round for
            checking secondary structure matches to make this easier.
*/
ZONE *ReadSSAP(FILE *fp)
{
   ZONE *zones = NULL,
        *z;
   int  resnum1, resnum2,
        score,
        start1 = 0, start2 = 0,
        end1 = 0,   end2 = 0;
   char buffer[MAXBUFF],
        aa1, aa2, str1, str2;
   

   while(fgets(buffer,MAXBUFF,fp))
   {
      TERMINATE(buffer);
      if(strlen_nospace(buffer))
      {
/* For older SSAP with no insert codes
         fsscanf(buffer,"%3d%1x%c%1x%c%2x%3d%2x%c%1x%c%1x%3d",
                 &resnum1, &str1, &aa1,
                 &score,
                 &aa2, &str2, &resnum2);
*/
         fsscanf(buffer,"%3d%1x%c%3x%c%2x%3d%2x%c%3x%c%1x%3d",
                 &resnum1, &str1, &aa1,
                 &score,
                 &aa2, &str2, &resnum2);


         /* If neither residue is an insert and both are E or both are H 
            or the -n flag has been set and neither are E or H then we
            are in a zone so record this fact
         */
         if((aa1 != ' ') && (aa2 != ' ') &&
            (((str1 == 'E') && (str2 == 'E')) ||
             ((str1 == 'H') && (str2 == 'H')) ||
             (gDoRandomCoil &&
              (str1 != 'E') && (str2 != 'E') &&
              (str1 != 'H') && (str2 != 'H'))))
         {
            if(!start1 || !start2)
            {
               start1 = resnum1;
               start2 = resnum2;
            }
            end1 = resnum1;
            end2 = resnum2;
         }
         else /* We've come out of a zone; store the last one           */
         {
            if(start1 && start2)
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

               z->start[0] = start1;
               z->start[1] = start2;
               z->end[0]   = end1;
               z->end[1]   = end2;

               start1 = start2 = 0;
            }
         }
      }
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
BOOL DefineCore(FILE *outfp, PDB *pdb1, PDB *pdb2, ZONE *zones, REAL dcut)
{
   int  count = 0,
        last  = 0,
        iter  = 0,
        natom1,
        natom2;
   REAL rm[3][3];
   PDB  *pdbca1,
        *pdbca2,
        **idx1,
        **idx2;
   
   /* Duplicate the PDB linked lists                                    */
   if((pdbca1 = blDupePDB(pdb1)) == NULL)
      return(FALSE);
   if((pdbca2 = blDupePDB(pdb2)) == NULL)
   {
      FREELIST(pdbca1,PDB);
      return(FALSE);
   }
   
   /* Reduce to CA only                                                 */
   pdbca1 = blSelectCaPDB(pdbca1);
   pdbca2 = blSelectCaPDB(pdbca2);
   
   SetBValByZone(pdbca1,zones,0);
   SetBValByZone(pdbca2,zones,1);

   count = CountCore(pdbca1);

   if((idx1 = blIndexPDB(pdbca1, &natom1))==NULL)
   {
      FREELIST(pdbca1,PDB);
      FREELIST(pdbca2,PDB);
      return(FALSE);
   }
   if((idx2 = blIndexPDB(pdbca2, &natom2))==NULL)
   {
      FREELIST(pdbca1,PDB);
      FREELIST(pdbca2,PDB);
      free(idx1);
      return(FALSE);
   }

   if(gInitialCut)
   {
      FitCaPDBBFlag(pdbca1, pdbca2, rm);
      if(!DoCut(idx1, natom1, idx2, natom2, zones, dcut*dcut))
         return(FALSE);
      count = CountCore(pdbca1);

      if(gVerbose)
      {
         fprintf(outfp,"\nCore after removing residues > 3.0A:\n");
         WriteTextOutput(outfp, zones);
      }
   }
   
   iter=0;
   while(last != count)
   {
      FitCaPDBBFlag(pdbca1, pdbca2, rm);
      UpdateBValues(idx1, natom1, idx2, natom2, zones, dcut*dcut);
      last = count;
      count = CountCore(pdbca1);
      if(++iter > MAXITER)
      {
         fprintf(stderr,"Warning: Maximum number of iterations (%d) \
exceeded!\n",MAXITER);
         break;
      }
   }
   
   free(idx1);
   free(idx2);

   return(TRUE);
}

/************************************************************************/
/*>BOOL DoCut(PDB **idx1, int natom1, PDB **idx2, int natom2,
              ZONE *zones, REAL cutsq)
   ----------------------------------------------------------
   Performs the initial cut of pairs which deviate by >3.0A

   06.12.96 Original   By: ACRM
*/
BOOL DoCut(PDB **idx1, int natom1, PDB **idx2, int natom2,
           ZONE *zones, REAL cutsq)
{
   ZONE *z, *zend, *znext;
   int  i, j,
        start1, end1,
        start2, end2;
   BOOL split,
        ok;
   
   for(z=zones; z!=NULL; NEXT(z))
   {
      /* Look for the start of the zone                                 */
      for(start1=0; start1<natom1; start1++)
      {
         if(idx1[start1]->resnum == z->start[0])
            break;
      }
      for(start2=0; start2<natom2; start2++)
      {
         if(idx2[start2]->resnum == z->start[1])
            break;
      }

      /* Loop for the end of the zone                                   */
      for(end1=0; end1<natom1; end1++)
      {
         if(idx1[end1]->resnum == z->end[0])
            break;
      }
      for(end2=0; end2<natom2; end2++)
      {
         if(idx2[end2]->resnum == z->end[1])
            break;
      }

      
      /* Step through the zone to see if we are within the cutoff       */
      split = FALSE;
      for(i=start1, j=start2; i<=end1 && j<=end2; i++, j++)
      {
         if(DISTSQ(idx1[i], idx2[j]) > cutsq)
         {
            split = TRUE;
            idx1[i]->bval = idx2[j]->bval = (REAL)0.0;
         }
      }

      /* See if there were any bits out of range which need the zones
         to be split
      */
      if(split)
      {
         /* See if the whole zone was out of range                      */
         ok = FALSE;
         for(i=start1, j=start2; i<=end1 && j<=end2; i++, j++)
         {
            if((idx1[i]->bval == (REAL)10.0) ||
               (idx2[j]->bval == (REAL)10.0))
            {
               ok = TRUE;
               break;
            }
         }
         /* If the zone contains no pairs within 3.0A, mark it for
            deletion by setting all values to -9999
         */
         if(!ok)
         {
            z->start[0] = z->start[1] = 
               z->end[0] = z->end[1] = (-9999);
         }
         else
         {
            /* There were some parts which are still required, but it's
               been modified, so we need to split the zones up
         
               First see if we've lost residues from the start of the
               zone
            */
            while((idx1[start1]->bval == (REAL)0.0) ||
                  (idx2[start2]->bval == (REAL)0.0))
            {
               start1++;
               start2++;
               z->start[0] = idx1[start1]->resnum;
               z->start[1] = idx2[start2]->resnum;
            }
            /* Now remove residues from the end of the zone in the same
               way
            */
            while((idx1[end1]->bval == (REAL)0.0) ||
                  (idx2[end2]->bval == (REAL)0.0))
            {
               end1--;
               end2--;
               z->end[0] = idx1[end1]->resnum;
               z->end[1] = idx2[end2]->resnum;
            }

            /* See if the new zone is split                             */
            split = FALSE;
            for(i=start1, j=start2; i<=end1 && j<=end2; i++, j++)
            {
               if((idx1[i]->bval == (REAL)0.0) ||
                  (idx2[j]->bval == (REAL)0.0))
               {
                  split = TRUE;
                  break;
               }
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
               
               zend->end[0] = z->end[0];
               zend->end[1] = z->end[1];
               /* Step back through the zone to find the start of this
                  subzone
               */
               while((idx1[end1]->bval == (REAL)10.0) &&
                     (idx2[end2]->bval == (REAL)10.0))
               {
                  end1--;
                  end2--;
                  zend->start[0] = idx1[end1]->resnum;
                  zend->start[1] = idx2[end2]->resnum;
               }
               /* Now step back to the end of the previous subzone      */
               while((idx1[end1]->bval == (REAL)0.0) ||
                     (idx2[end2]->bval == (REAL)0.0))
               {
                  end1--;
                  end2--;
                  z->end[0] = idx1[end1]->resnum;
                  z->end[1] = idx2[end2]->resnum;
               }

               /* Test again to see if it's split                       */
               split = FALSE;
               for(i=start1, j=start2; i<=end1 && j<=end2; i++, j++)
               {
                  if((idx1[i]->bval == (REAL)0.0) ||
                     (idx2[j]->bval == (REAL)0.0))
                  {
                     split = TRUE;
                     break;
                  }
               }
            }
         }
      }
   }
      
   return(TRUE);
}

/************************************************************************/
/*>void UpdateBValues(PDB **idx1, int natom1, PDB **idx2, int natom2,
                      ZONE *zones, REAL cutsq)
   ------------------------------------------------------------------
   Update the B-values and the current zones by extending out from
   the secondary structure regions

   14.11.96 Original   By: ACRM
   13.03.97 Now checks that residues are not already in a zone before
            adding them to the current zone. Fixes a problem at the
            zone-merge stage where the merged zones could end up with
            different numbers of residues.
*/
void UpdateBValues(PDB **idx1, int natom1, PDB **idx2, int natom2,
                   ZONE *zones, REAL cutsq)
{
   ZONE *z;
   int  i, j;

   
   for(z=zones; z!=NULL; NEXT(z))
   {
      if(z->start[0] < -9998)
         continue;
      
      /* Look for the start of the zone                                 */
      for(i=0; i<natom1; i++)
      {
         if(idx1[i]->resnum == z->start[0])
            break;
      }
      for(j=0; j<natom2; j++)
      {
         if(idx2[j]->resnum == z->start[1])
            break;
      }
      
      /* Step back from the start seeing if we are within the cutoff    */
      i--; j--;
      while(i>=0 && j>=0)
      {
         if((DISTSQ(idx1[i], idx2[j]) > cutsq) ||
            (idx1[i]->bval > (REAL)5.0)        ||
            (idx2[j]->bval > (REAL)5.0))
            break;
         else
            idx1[i]->bval = idx2[j]->bval = (REAL)10.0;
            
         i--; j--;
      }
      i++; j++;

      z->start[0] = idx1[i]->resnum;
      z->start[1] = idx2[j]->resnum;
      
      /* Look for the end of the zone                                   */
      for(i=0; i<natom1; i++)
      {
         if(idx1[i]->resnum == z->end[0])
            break;
      }
      for(j=0; j<natom2; j++)
      {
         if(idx2[j]->resnum == z->end[1])
            break;
      }
      
      /* Step forward from the end seeing if we are within the cutoff   */
      i++; j++;
      while(i<natom1 && j<natom2)
      {
         if((DISTSQ(idx1[i], idx2[j]) > cutsq) ||
            (idx1[i]->bval > (REAL)5.0)        ||
            (idx2[j]->bval > (REAL)5.0))
            break;
         else
            idx1[i]->bval = idx2[j]->bval = (REAL)10.0;
         i++; j++;
      }
      i--; j--;
      
      z->end[0] = idx1[i]->resnum;
      z->end[1] = idx2[j]->resnum;
   }
}

/************************************************************************/
/*>void SetBValByZone(PDB *pdb, ZONE *zones, int which)
   ----------------------------------------------------
   Set B-values to 10 if in the zones, otherwise to 0.0

   14.11.96 Original   By: ACRM
*/
void SetBValByZone(PDB *pdb, ZONE *zones, int which)
{
   PDB  *p;
   ZONE *z;

   for(p=pdb; p!=NULL; NEXT(p))
      p->bval = (REAL)0.0;
   
   for(p=pdb; p!=NULL; NEXT(p))
   {
      for(z=zones; z!=NULL; NEXT(z))
      {
         if((p->resnum >= z->start[which]) &&
            (p->resnum <= z->end[which]))
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




/************************************************************************/
/*>void WriteTextOutput(FILE *fp, ZONE *zones)
   -------------------------------------------
   Writes the zones out in text format

   14.11.96 Original   By: ACRM
   06.12.96 Added check that zones have not been blanked out
*/
void WriteTextOutput(FILE *fp, ZONE *zones)
{
   ZONE *z;
   
   for(z=zones; z!=NULL; NEXT(z))
   {
      if(z->start[0] > -9999)
      {
         fprintf(fp,"%d-%d : %d-%d\n",
                 z->start[0],z->end[0],
                 z->start[1],z->end[1]);
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
   26.06.02 V1.4
   05.11.25 V1.5
*/
void Usage(void)
{
   fprintf(stderr,"\nFindCore V1.5 (c) 1996-2025, Prof. Andrew C.R. Martin, \
UCL.\n");

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



   fprintf(stderr,"\nFindCore defines a protein structurally conserved \
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
/*>ZONE *MergeZones(ZONE *zones)
   -----------------------------
   Merges zones in the zone linked list of they overlap and removes zones
   which are subsets of other zones

   14.11.96 Original   By: ACRM
   06.12.96 Removes zones marked for deletion
   13.03.97 Added merging of abutting zones (required since change to
            zone creation where residues not added to a zone if they
            are already in another zone stops them from being subsets)
   26.06.02 Fixed bug in deleting zones
*/
ZONE *MergeZones(ZONE *zones)
{
   ZONE *z, *zt;
   BOOL finished = FALSE;

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
         if((z->end[0] >= z->next->start[0]) ||
            (z->end[1] >= z->next->start[1]))
         {
            if((z->end[0] != z->next->end[0]) &&
               (z->end[1] != z->next->end[1]))
            {
               finished = FALSE;
               z->end[0] = z->next->end[0];
               z->end[1] = z->next->end[1];
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
         /* If these two zones are both abutting                        */
         if((z->end[0]+1 == z->next->start[0]) &&
            (z->end[1]+1 == z->next->start[1]))
         {
            /* Reset start of second zone to start of first zone        */
            z->next->start[0] = z->start[0];
            z->next->start[1] = z->start[1];
            
            /* Remove z from the linked list                            */
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


