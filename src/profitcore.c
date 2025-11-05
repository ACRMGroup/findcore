#include <stdio.h>
#include <stdlib.h>
#include "bioplib/macros.h"
#include "bioplib/pdb.h"

/************************************************************************/
typedef struct _zone
{
   int  start[2],
        stop[2];
   char startresid[2][16],
        stopresid[2][16];
   struct _zone *next;
} ZONE;

#define MAXBUFF 320

/************************************************************************/
ZONE *ReadProFitZones(FILE *fp);
void PrintZones(FILE *fp, ZONE *zones);
BOOL MapZones(ZONE *zones, int strucNum, PDB *pdb);
void AnnotateZones(ZONE *zones, int strucNum, PDB *pdb);
void Die(char *msg, char *submsg, int status);

/************************************************************************/
int main(int argc, char **argv)
{
   ZONE *zones = NULL;
   FILE *fp, *fpP1, *fpP2;
   int natoms;
   PDB *pdb1, *pdb2;

   if((fp = fopen("zones.txt", "r"))==NULL)
      Die("Unable to open zones file: ", "zones.txt", 1);

   if((fpP1 = fopen("pdb1yqv_0P.mar", "r"))==NULL)
      Die("Unable to open first PDB input file: ", "pdb1yqv_0P.mar", 1);

   if((pdb1 = blReadPDB(fpP1, &natoms))==NULL)
      Die("No atoms read from first PDB input file: ", "pdb1yqv_0P.mar", 1);

   if((fpP2 = fopen("pdb8fab_0.mar", "r"))==NULL)
      Die("Unable to open second PDB input file: ", "pdb8fab_0.mar", 1);
      
   if((pdb2 = blReadPDB(fpP2, &natoms))==NULL)
      Die("No atoms read from second PDB input file: ", "pdb8fab_0.mar", 1);
   
   if((zones = ReadProFitZones(fp))==NULL)
      Die("Unable to read zones from the zones file", NULL, 1);
   
   if(!MapZones(zones, 0, pdb1))
      Die("No memory for mapping zones", NULL, 1);

   if(!MapZones(zones, 1, pdb2))
      Die("No memory for mapping zones", NULL, 1);
      
   AnnotateZones(zones, 0, pdb1);
   AnnotateZones(zones, 1, pdb2);
   
   PrintZones(stdout, zones);
   
   return(0);
}


/************************************************************************/
void PrintZones(FILE *fp, ZONE *zones)
{
   ZONE *z;
   
   for(z=zones; z!=NULL; NEXT(z))
   {
      fprintf(fp, "%s (%d) to %s (%d) with %s (%d) to %s (%d)\n",
              z->startresid[0], z->start[0],
              z->stopresid[0],  z->stop[0],
              z->startresid[1], z->start[1],
              z->stopresid[1],  z->stop[1]);
   }
}


/************************************************************************/
ZONE *ReadProFitZones(FILE *fp)
{
   ZONE *zones = NULL,
        *z     = NULL;
   char buffer[MAXBUFF],
        junk[16];
   int  start1, stop1, start2, stop2;
   
   while(fgets(buffer, MAXBUFF, fp))
   {
      if(sscanf(buffer, "%d %s %d %s %d %s %d", 
                &start1, junk,
                &stop1,  junk,
                &start2, junk,
                &stop2) == 7)
      {
         /* Allocate next position in linked list */
         if(zones == NULL)
         {
            INIT(zones, ZONE);
            z=zones;
         }
         else
         {
            ALLOCNEXT(z, ZONE);
         }
         if(z==NULL)
         {
            FREELIST(zones, ZONE);
            return(NULL);
         }
         
         /* Populate linked list */
         z->start[0] = start1;
         z->stop[0]  = stop1;
         z->start[1] = start2;
         z->stop[1]  = stop2;
      }
   }
   return(zones);
}


/************************************************************************/
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
}


/************************************************************************/
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
void Die(char *msg, char *submsg, int status)
{
   if(submsg != NULL)
      fprintf(stderr, "Error (zones2ssap) %s%s\n",msg, submsg);
   else
      fprintf(stderr, "Error (zones2ssap) %s\n",msg);

   exit(status);
}
