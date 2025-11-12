/************************************************************************/
/**

   Program:    profitcore
   \file       profitcore.c
   
   \version    V1.1
   \date       12.11.25   
   \brief      Identify protein core from ProFit iterative fit
   
   \copyright  (c) Prof Andrew C. R. Martin 2025
   \author     Prof. Andrew C. R. Martin
   \par
               abYinformatics, Ltd
               www.bioinf.org.uk
   \par
               andrew@bioinf.org.uk
               andrew@abyinformatics.com
               
**************************************************************************

   This code is released under the GPL V3.0

**************************************************************************

   Description:
   ============

**************************************************************************

   Usage:
   ======

**************************************************************************

   Revision History:
   =================
-  V1.0   05.11.25  Original   By: ACRM
-  V1.1   12.11.25  Added multiple fitting (-m)

*************************************************************************/
/* Includes
*/
#include <stdio.h>
#include <stdlib.h>
#include "bioplib/macros.h"
#include "bioplib/pdb.h"

/************************************************************************/
/* Defines and macros
*/
typedef struct _zone
{
   int  start[2],
        stop[2];
   char startresid[2][16],
        stopresid[2][16];
   struct _zone *next;
} ZONE;
typedef struct _mzone
{
   ZONE *zones;
   struct _mzone *next;
}  MZONE;

#define MAXBUFF 320

/************************************************************************/
/* Globals
*/

/************************************************************************/
/* Prototypes
*/
MZONE *ReadProFitZones(FILE *fp);
void PrintZones(FILE *fp, ZONE *zones);
BOOL MapZones(ZONE *zones, int strucNum, PDB *pdb);
void AnnotateZones(ZONE *zones, int strucNum, PDB *pdb);
BOOL WriteFile(PDB *pdb, char *filename);
void Die(char *msg, char *submsg, int status);
void Usage(void);
BOOL ParseCmdLine(int argc, char **argv, char *zoneFile,
                  char *pdbFile1, char *pdbFile2,
                  char *outFile1, char *outFile2,
                  BOOL *multi,    char *multiFile,
                  char *extIn,    char *extOut);

/************************************************************************/
/*>int main(int argc, char **argv)
   -------------------------------
*//**
   Main program for core finding
   
-  05.11.25 Original   By: ACRM
*/
int main(int argc, char **argv)
{
   MZONE *mzones = NULL;
   FILE  *fp, *fpP1, *fpP2;
   int   natoms;
   PDB   *pdb1, *pdb2;
   char  zoneFile[MAXBUFF],
         pdbFile1[MAXBUFF],
         pdbFile2[MAXBUFF],
         outFile1[MAXBUFF],
         outFile2[MAXBUFF],
         multiFile[MAXBUFF],
         extIn[16],
         extOut[16];
   BOOL  multi = FALSE;
   
   if(ParseCmdLine(argc, argv, zoneFile, pdbFile1, pdbFile2,
                   outFile1, outFile2, &multi, multiFile,
                   extIn, extOut))
   {
      if((fp = fopen(zoneFile, "r"))==NULL)
         Die("Unable to open zones file: ", zoneFile, 1);
      if((mzones = ReadProFitZones(fp))==NULL)
         Die("Unable to read zones from the zones file", NULL, 1);

      if(multi)
      {
         fprintf(stderr, "\n\n**** WRITE ME ****\n\n");
      }
      else
      {
         if((fpP1 = fopen(pdbFile1, "r"))==NULL)
            Die("Unable to open first PDB input file: ",
                pdbFile1, 1);
         if((pdb1 = blReadPDB(fpP1, &natoms))==NULL)
            Die("No atoms read from first PDB input file: ",
                pdbFile1, 1);
         
         if((fpP2 = fopen(pdbFile2, "r"))==NULL)
            Die("Unable to open second PDB input file: ",
                pdbFile2, 1);
         if((pdb2 = blReadPDB(fpP2, &natoms))==NULL)
            Die("No atoms read from second PDB input file: ",
                pdbFile2, 1);
         if(!MapZones(mzones->zones, 0, pdb1))
            Die("No memory for mapping zones", NULL, 1);
         if(!MapZones(mzones->zones, 1, pdb2))
            Die("No memory for mapping zones", NULL, 1);

         PrintZones(stdout, mzones->zones);

         if(outFile1[0] != '\0')
         {
            AnnotateZones(mzones->zones, 0, pdb1);
            if(!WriteFile(pdb1, outFile1))
               Die("Unable to write PDB file: ", outFile1, 1);
         }
         
         if(outFile2[0] != '\0')
         {
            AnnotateZones(mzones->zones, 1, pdb2);
            if(!WriteFile(pdb2, outFile2))
               Die("Unable to write PDB file: ", outFile2, 1);
         }
      }
   }
   else
   {
      Usage();
   }
   
   return(0);
}


/************************************************************************/
/*>void PrintZones(FILE *fp, ZONE *zones)
   --------------------------------------
*//**

   \param[in]     *fp    File pointer for printing
   \param[in]     *zones Linked list of zones

   Prints the converted zone information
   
-  05.11.25 Original   By: ACRM
*/
void PrintZones(FILE *fp, ZONE *zones)
{
   ZONE *z;
   
   for(z=zones; z!=NULL; NEXT(z))
   {
      fprintf(fp, "%s to %s with %s to %s\n",
              z->startresid[0],
              z->stopresid[0], 
              z->startresid[1],
              z->stopresid[1]);
   }
}


/************************************************************************/
/*>MZONE *ReadProFitZones(FILE *fp)
   -------------------------------
*//**

   \param[in]    *fp    File pointer for zones file

   Reads the zone information from the file which is simply cut and
   paste from the ProFit status message

-  05.11.25 Original   By: ACRM
-  12.11.25 Now returns a list of multiple zones
            i.e. a linked list of linked lists zones
*/
MZONE *ReadProFitZones(FILE *fp)
{
   MZONE *mzones = NULL,
         *mz     = NULL;
   ZONE  *z     = NULL;
   char  buffer[MAXBUFF],
         junk[16];
   int   start1, stop1, start2, stop2;

   INIT(mzones, MZONE);
   mz=mzones;
   if(mz==NULL)
      return(NULL);
   
   while(fgets(buffer, MAXBUFF, fp))
   {
      if(sscanf(buffer, "%d %s %d %s %d %s %d", 
                &start1, junk,
                &stop1,  junk,
                &start2, junk,
                &stop2) == 7)
      {
         /* Allocate next position in linked list */
         if(mz->zones == NULL)
         {
            INIT(mz->zones, ZONE);
            z=mz->zones;
         }
         else
         {
            ALLOCNEXT(z, ZONE);
         }
         if(z==NULL)
         {
            FREELIST(mz->zones, ZONE);
            return(NULL);
         }
         
         /* Populate linked list */
         z->start[0] = start1;
         z->stop[0]  = stop1;
         z->start[1] = start2;
         z->stop[1]  = stop2;
      }
   }
   return(mzones);
}


/************************************************************************/
/*>BOOL MapZones(ZONE *zones, int strucNum, PDB *pdb)
   --------------------------------------------------
*//**

   \param[in,out] *zones    Linked list of zones
   \param[in]     strucNum  The structure being mapped (0 or 1)
   \param[in]     pdb       PDB linked list for the structure

   Maps the sequentially numbered zones to PDB residue IDs
   
-  05.11.25 Original   By: ACRM
*/
BOOL MapZones(ZONE *zones, int strucNum, PDB *pdb)
{
   ZONE *z;
   PDB  **idx  = NULL,
        *pdbca = NULL;
   int  natoms;
   char *sel[2];
   SELECT(sel[0], "CA  ");

   /* Create new PDB list of only C-alphas                              */
   if((pdbca = blSelectAtomsPDBAsCopy(pdb, 1, sel, &natoms))==NULL)
      return(FALSE);

   /* Index the PDB linked list                                         */
   if((idx = blIndexPDB(pdbca, &natoms))==NULL)
   {
      FREELIST(pdbca, PDB);
      return(FALSE);
   }

   /* Use the sequential residue counts in the CA-only list to identify
      the residues and put them in our zone list
   */
   for(z=zones; z!=NULL; NEXT(z))
   {
      PDB *p = NULL;

      p = idx[z->start[strucNum] - 1];
      MAKERESID(z->startresid[strucNum], p);

      p = idx[z->stop[strucNum] - 1];
      MAKERESID(z->stopresid[strucNum], p);
   }

   FREE(idx);
   FREELIST(pdbca, PDB);

   return(TRUE);
}


/************************************************************************/
/*>void AnnotateZones(ZONE *zones, int strucNum, PDB *pdb)
   -------------------------------------------------------
*//**

   \param[in]     *zones    Linked list of zones after mapping
   \param[in]     strucNum  The structure being mapped (0 or 1)
   \param[in,out] *pdb      PDB linked list     

   Updates the temperature factor column in the PDB file such that
   all atoms are initially set to zero and then those in zones are
   set to 1.0

-  05.11.25 Original   By: ACRM
*/
void AnnotateZones(ZONE *zones, int strucNum, PDB *pdb)
{
   ZONE *z;
   PDB  *p;

   /* Set all B-values to zero                                          */
   for(p=pdb; p!=NULL; NEXT(p))
      p->bval = 0;

   /* Find the zone residues in the PDB list and set the B-values to
      indicate them as of interest
   */
   for(z=zones; z!=NULL; NEXT(z))
   {
      PDB *pStart = NULL,
          *pStop  = NULL;

      /* Find this residue in the PDB linked list                       */
      pStart = blFindResidueSpec(pdb, (char *)z->startresid);
      pStop  = blFindResidueSpec(pdb, (char *)z->stopresid);
      pStop  = blFindNextResidue(pStop);

      /* Set all atoms to have B-value of one                           */
      for(p=pStart; p!=pStop; NEXT(p))
         p->bval = 1;
   }
}


/************************************************************************/
/*>BOOL WriteFile(PDB *pdb, char *filename)
   ----------------------------------------
*//**

   \param[in]     *pdb       PDB linked list
   \param[in]     *filename  Name of output file to create
   
   Simple wrapper to open a file and write a PDB linked list to it
   
-  05.11.25 Original   By: ACRM
*/
BOOL WriteFile(PDB *pdb, char *filename)
{
   FILE *fp = NULL;
   
   if((fp=fopen(filename, "w"))!=NULL)
   {
      blWritePDB(fp, pdb);
      fclose(fp);
      return(TRUE);
   }
   return(FALSE);
}


/************************************************************************/
/*>void Die(char *msg, char *submsg, int status)
   ---------------------------------------------
*//**

   \param[in]     *msg     Main error message
   \param[in]     *submsg  Optional subsiduary message
   \param[in]     status   Exit status
   
   Die with error message

-  05.11.25 Original   By: ACRM
*/
void Die(char *msg, char *submsg, int status)
{
   if(submsg != NULL)
      fprintf(stderr, "Error (zones2ssap) %s%s\n",msg, submsg);
   else
      fprintf(stderr, "Error (zones2ssap) %s\n",msg);

   exit(status);
}

/************************************************************************/
/*>void Usage(void)
   ----------------
*//**

   Usage message

-  05.11.25 Original   By: ACRM
-  12.11.25 Added -m
*/
void Usage(void)
{
   printf("\nprofitcore V1.1 (c) 2025, Prof Andrew C.R. Martin, \
abYinformatics\n");
   printf("\nUsage: profitcore [-o1 file] [-o2 file] zoneFile \
pdbfile1 pdbfile2\n");
   printf("  -or-\n");    
   printf("       profitcore [-m multiFile] [-xi inExt] [-xo outExt] \
zoneFile \n");
   printf("\n");

   printf("       -o1 Specify first output PDB file\n");
   printf("       -o2 Specify second output PDB file\n");
   printf("\nMulti-mode\n");
   printf("       -m  Specify the multi file as used by ProFit\n");
   printf("       -xi Specify input extension for fitted files\n");
   printf("       -xo Specify output extension for core PDB files\n");

   printf("profitcore converts the sequentially numbered zones output \
by ProFit into\n");
   printf("residue identifiers and, optionally, generates PDB files \
with the B-value\n");
   printf("used to indicate residues that are in those zones.\n");
   printf("\n");
   printf("By running ProFit with the commands:\n");
   printf("   ALIGN\n");
   printf("   ITER\n");
   printf("   FIT\n");
   printf("   STATUS\n");
   printf("it will perform an iterative structural alignment creating \
zones based\n");
   printf("on a dynamic programming alignment discarding residue pairs \
with C-alphas\n");
   printf("more than 3A apart, thus identifying a conserved core. You \
can alter the\n");
   printf("3A threshold by providing a distance as a parameter to the \
ITER command.\n");
   printf("Thus ITER 2.0 would identify a stricter core, while ITER 4.0 \
would\n");
   printf("allow more flexibility.\n\n");
}

/************************************************************************/
/*>BOOL ParseCmdLine(int argc, char **argv, char *zoneFile,
                  char *pdbFile1, char *pdbFile2,
                  char *outFile1, char *outFile2,
                  BOOL *multi,    char *multiFile,
                  char *extIn,    char *extOut)
   --------------------------------------------------------
*//**

   \param[in]     argc      Argument count
   \param[in]     argv      Arguments
   \param[out]    zoneFile  The name of the zone file
   \param[out]    pdbFile1  The name of the 1st PDB file
   \param[out]    pdbFile2  The name of the 2nd PDB file
   \param[out]    outFile1  The name of the (optional) 1st output file
   \param[out]    outFile2  The name of the (optional) 2nd output file
   \param[out]    multi     We are finding the core from multiple files
   \param[out]    multiFile The name of the multi-fitting file
   \param[out]    extIn     The extension of the multiple fitted PDB files
   \param[out]    extOut    The extension of the multiple core PDB files

   Parses the command line

-  05.11.25 Original   By: ACRM
-  12.11.25 Added -m and -x
*/
BOOL ParseCmdLine(int argc, char **argv, char *zoneFile,
                  char *pdbFile1, char *pdbFile2,
                  char *outFile1, char *outFile2,
                  BOOL *multi,    char *multiFile,
                  char *extIn,    char *extOut)
{
   BOOL gotFile  = FALSE;
   argc--;
   argv++;

   zoneFile[0]     =
      pdbFile1[0]  = pdbFile2[0] =
      outFile1[0]  = outFile2[0] =
      multiFile[0] = '\0';
   strcpy(extIn,  "fit");
   strcpy(extOut, "cor");
   *multi = FALSE;

   while(argc)
   {
      if(argv[0][0] == '-')
      {
         switch(argv[0][1])
         {
         case 'o':
            if(*multi)
            {
               fprintf(stderr,"\nError! You cannot mix -o with \
-x or -m\n\n");
               return(FALSE);
            }
            gotFile = TRUE;
            
            switch(argv[0][2])
            {
            case '1':
               argc--; argv++;
               strcpy(outFile1, argv[0]);
               break;
            case '2':
               argc--; argv++;
               strcpy(outFile2, argv[0]);
               break;
            default:
               return(FALSE);
            }
            break;
         case 'x':
            if(gotFile)
            {
               fprintf(stderr,"\nError! You cannot mix -o with -x \
or -m\n\n");
               return(FALSE);
            }
            *multi = TRUE;
            
            switch(argv[0][2])
            {
            case 'i':
               argc--; argv++;
               strcpy(extIn, argv[0]);
               break;
            case 'o':
               argc--; argv++;
               strcpy(extOut, argv[0]);
               break;
            default:
               return(FALSE);
            }
            break;
         case 'm':
            if(gotFile)
            {
               fprintf(stderr,"\nError! You cannot mix -o with -m\n\n");
               return(FALSE);
            }
            *multi = TRUE;
            argc--; argv++;
            strcpy(multiFile, argv[0]);
            break;
         case 'h':
         default:
            return(FALSE);
         }
      }
      else
      {
         if(*multi)
         {
            /* Check that there is 1 argument left                      */
            if(argc != 1)
               return(FALSE);
            /* Copy the file names                                      */
            strcpy(zoneFile, argv[0]);
         }
         else
         {
            /* Check that there are 3 arguments left                    */
            if(argc != 3)
               return(FALSE);
            /* Copy the file names                                      */
            strcpy(zoneFile, argv[0]);
            strcpy(pdbFile1, argv[1]);
            strcpy(pdbFile2, argv[2]);
         }
         return(TRUE);
      }
      argc--; argv++;
   }
   return(FALSE);
}
