/*----------------------------------------------------------------------*/
/*      eccdemo8.c      ecc demo program  - version 1.7                 */
/*                                                                      */
/*      (c) 2017 by Jeff Reid                                           */
/*----------------------------------------------------------------------*/
#define _CRT_SECURE_NO_WARNINGS 1       /* disable sscanf warnings */

#include <stdio.h>
#include <string.h>
#include <memory.h>

typedef unsigned  char  BYTE;

#define DISPLAYE 0                      /* display euclid stuff */
#define DISPLAYM 0                      /* display matrix stuff */
#define DISPLAYI 0                      /* display matrix inv */

/*                                         maximum # parity bytes */
#define MAXPAR 20

typedef struct{                         /* vector structure */
    BYTE  size;
    BYTE  data[255];
}VECTOR;

typedef struct{                         /* matrix structure */
    BYTE  nrows;
    BYTE  ncols;
    BYTE  data[(MAXPAR)*(MAXPAR)*2];
}MATRIX;

typedef struct{                         /* euclid structure */
    BYTE  size;                         /* # of data bytes */
    BYTE  indx;                         /* index to right side */
    BYTE  data[MAXPAR+2];               /* left and right side data */
}EUCLID;

/*----------------------------------------------------------------------*/
/*      data                                                            */
/*----------------------------------------------------------------------*/
static int aiGA[60] =                   /* GF's and Alpha's */
   {0x11b,0x03,0x11d,0x02,0x12b,0x02,0x12d,0x02,
    0x139,0x03,0x13f,0x03,0x14d,0x02,0x15f,0x02,
    0x163,0x02,0x165,0x02,0x169,0x02,0x171,0x02,
    0x177,0x03,0x17b,0x0b,0x187,0x02,0x18b,0x0b,
    0x18d,0x02,0x19f,0x03,0x1a3,0x03,0x1a9,0x02,
    0x1b1,0x07,0x1bd,0x07,0x1c3,0x02,0x1cf,0x02,
    0x1d7,0x07,0x1dd,0x07,0x1e7,0x02,0x1f3,0x0d,
    0x1f5,0x02,0x1f9,0x03};

static BYTE     abExp[512];             /* antilog table */
static BYTE     abLog[256];             /* log table */

static BYTE     abUser[80];             /* user input buffer */

/*  pErrors = poly wtih roots Alpha^location == pLambda reversed */
/*  pLambda = poly with roots 1/(Alpha^location) */
/*  some articles use the name simga instead of lambda */

static VECTOR   vGenRoots;              /* generator roots */
static VECTOR   pGenPoly;               /* generator poly */
static VECTOR   vData;                  /* data */
static VECTOR   vErsf;                  /* erasure flags */
static VECTOR   vParities;              /* parities */
static VECTOR   vSyndromes;             /* syndromes */
static VECTOR   vErsLocators;           /* erasure locators */
static VECTOR   pErasures;              /* erasure poly */
static VECTOR   vModSyndromes;          /* modified syndromes */
static MATRIX   mModSyndromes;          /* mod syndrome matrix */
static MATRIX   imModSyndromes;         /* inv mod syndrome matrix */
static VECTOR   pErrors;                /* error locator poly */
static VECTOR   vErrLocators;           /* error locators */
static VECTOR   vLocators;              /* erasure + error locators */
static VECTOR   pLocators;              /* locator poly */
static MATRIX   mLocators;              /* locator matrix */
static MATRIX   imLocators;             /* inverse locator matrix */
static MATRIX   mModLocators;           /* Matrix of modified Locators */
static VECTOR   vErrValues;             /* error values */
static VECTOR   vInvLocators;           /* inverse (1/x) locators */
static VECTOR   vLocations;             /* locations n-1 to 0 */
static VECTOR   vOffsets;               /* offsets 0 to n-1 */

static EUCLID   E0;                     /* used by GenpErrorsE */
static EUCLID   E1;

static VECTOR   vB;                     /* used by GenpErrorsB */
static VECTOR   vC;
static VECTOR   vT;
static VECTOR   vBx;

static VECTOR   pLambda;                /* error 1/locator poly */
static VECTOR   pDLambda;               /* "derivative" of pLambda */
static VECTOR   pOmega;                 /* error evaluator poly */
static VECTOR   vForney;                /* error values from Forney */

static VECTOR   vErrValues;

static MATRIX   mAltInvLoc;
static VECTOR   vAltErrValues;

static BYTE     abId[MAXPAR*2+1];       /* array for MatrixInv */

static int      iNumData;               /* number data - parities */
static int      iGF;                    /* Galios Field Polynomial */
static BYTE     fSelfRcp;               /* Self-reciprocal Poly flag */
static BYTE     bAlpha;                 /* Alpha for this field */
static BYTE     bFCR;                   /* first consecutive root */
static BYTE     bLogFCR;                /* log first concecutive root */
static BYTE     bForneyC;               /* Forney correction factor */

/*----------------------------------------------------------------------*/
/*      code                                                            */
/*----------------------------------------------------------------------*/
static void     Encode(void);
static void     Decode(void);
static void     GenSyndromes(void);
static int      GenErasures(void);
static void     GenModSyndromes(void);
static void     GenpErrorsM(void);
static void     GenpErrorsE(void);
static void     GenpErrorsB(void);
static void     MergeLocators(void);
static void     GenLambda(void);
static void     GenOffsets(void);
static void     GenmLocators(void);
static void     GenErrValues(void);
static void     GenAltInvLoc(void);
static void     GenOmega(void);
static void     GenForneyErr(void);
static void     FixvData(void);
static int      Poly2Root(VECTOR *, VECTOR *);
static void     Root2Poly(VECTOR *, VECTOR *);
static int      MatrixInv(MATRIX *, MATRIX *);
static BYTE     GFAdd(BYTE, BYTE);
static BYTE     GFSub(BYTE, BYTE);
static BYTE     GFMpy(BYTE, BYTE);
static BYTE     GFDiv(BYTE, BYTE);
static BYTE     GFPow(BYTE, BYTE);
static void     InitGF(void);
static void     ShowEuclid(EUCLID *);
static void     ShowVector(VECTOR *);
static void     ShowMatrix(MATRIX *);
static int      Conrs(void);
static void     DoUser(void);

/*----------------------------------------------------------------------*/
/*      Encode                                                          */
/*----------------------------------------------------------------------*/
static void Encode(void)
{
BYTE    bQuot;                          /* quotient byte */
BYTE    bRem0, bRem1;                   /* partial remainders */
int     i, j;

    vParities.size = vGenRoots.size;    /* init vParities */
    memset(vParities.data, 0, vParities.size);

    for(j = 0; j < iNumData; j++){      /* generate vParities */
        bQuot = GFAdd(vData.data[j], vParities.data[0]);
        bRem0 = 0;

        for(i = vGenRoots.size; i;){
            bRem1 = GFSub(bRem0, GFMpy(bQuot, pGenPoly.data[i]));
            i--;
            bRem0 = vParities.data[i];
            vParities.data[i] = bRem1;}}

    for(i = 0; i < vParities.size; i++){    /* append parities */
        vData.data[iNumData+i] = GFSub(0, vParities.data[i]);}
}

/*----------------------------------------------------------------------*/
/*      Decode                                                          */
/*----------------------------------------------------------------------*/
static void Decode(void)
{
    GenSyndromes();             /* generate syndromes */
    if(GenErasures())           /* generate erasure stuff */
        return;
    GenModSyndromes();          /* generate modified syndromes */
    GenpErrorsM();              /* generate error poly */
    GenpErrorsE();              /* generate error poly */
    GenpErrorsB();              /* generate error poly */
    printf("pLambda (B):    ");
    ShowVector(&vC);
    MergeLocators();            /* merge erasures + errors */
    printf("pLocators:      ");
    ShowVector(&pLocators);
    GenLambda();                /* generate lambda */
    printf("pLambda:        ");
    ShowVector(&pLambda);
    printf("pDLambda:       ");
    ShowVector(&pDLambda);
    GenOffsets();               /* generate offsets */
    GenmLocators();             /* generate locator matrix */
    MatrixInv(&imLocators, &mLocators); /* invert */
    GenErrValues();             /* generate error values */
    GenAltInvLoc();             /* generate alt matrix, err values */
    GenOmega();                 /* generate Omega */
    printf("pOmega:         ");
    ShowVector(&pOmega);
    GenForneyErr();             /* generate forney err values */

    printf("vLocators:      ");
    ShowVector(&vLocators);
    printf("vInvLocators:   ");
    ShowVector(&vInvLocators);
    printf("imLocators:\n");
    ShowMatrix(&imLocators);
    printf("mAltInvLoc:\n");
    ShowMatrix(&mAltInvLoc);
    printf("vErrValues:     ");
    ShowVector(&vErrValues);
    printf("vAltErrValues:  ");
    ShowVector(&vAltErrValues);
    printf("vForney:        ");
    ShowVector(&vForney);
    printf("vOffsets:       ");
    ShowVector(&vOffsets);
    FixvData();
}

/*----------------------------------------------------------------------*/
/*      GenSyndromes    generate standard RS syndromes                  */
/*----------------------------------------------------------------------*/
static void GenSyndromes(void)
{
BYTE i,j;

    vSyndromes.size = vGenRoots.size;   /* set size */

    for(j = 0; j < vGenRoots.size; j++){
        vSyndromes.data[j] = vData.data[0]; /* generate a syndrome */
        for(i = 1; i < vData.size;  i++){
            vSyndromes.data[j] = GFAdd(vData.data[i],
                GFMpy(vGenRoots.data[j], vSyndromes.data[j]));}}

    printf("vSyndromes:     ");
    ShowVector(&vSyndromes);
}

/*----------------------------------------------------------------------*/
/*      GenErasures     generate vErsLocat...and pErasures              */
/*                                                                      */
/*      Scan vErsf right to left; for each non-zero flag byte,          */
/*      set vErsLocators to Alpha**power of corresponding location.     */
/*      Then convert these locators into a polynomial.                  */
/*----------------------------------------------------------------------*/
static int GenErasures(void)
{
int     i, j;
BYTE    bLcr;                           /* locator */

    j = 0;                              /* j = index into VErs... */
    bLcr = 1;                           /* init locator */
    for(i = vErsf.size; i;){
        i--;
        if(vErsf.data[i]){              /* if byte flagged */
            vErsLocators.data[j] = bLcr; /*     set locator */
            j++;
            if(j > vGenRoots.size){      /*    exit if too many */
                return(1);}}
        bLcr = GFMpy(bLcr, bAlpha);}     /* bump locator */

    vErsLocators.size = j;              /* set size */

    Root2Poly(&pErasures, &vErsLocators); /* generate poly */

    return(0);
}

/*----------------------------------------------------------------------*/
/*      GenModSyndromes         Generate VModSyndromes                  */
/*                                                                      */
/*      First step in finding unknown locations.                        */
/*      Reset all unknown location parameters.                          */
/*      Each modified syndrome = sum(syndromes[i++]*perasures[j--]      */
/*----------------------------------------------------------------------*/
static void     GenModSyndromes(void)
{
int     iM;                             /* vModSyn.. index */
int     iS;                             /* vSyn.. index */
int     iP;                             /* pErs.. index */

    vModSyndromes.size = vSyndromes.size-vErsLocators.size;

    if(vErsLocators.size == 0){         /* if no erasures, copy */
        memcpy(vModSyndromes.data, vSyndromes.data, 
               vSyndromes.size*sizeof(BYTE));
        goto genms0;}

    iS = 0;
    for(iM = 0; iM < vModSyndromes.size; iM++){ /* modify */
        vModSyndromes.data[iM] = 0;
        iP = pErasures.size;
        while(iP){
            iP--;
            vModSyndromes.data[iM] = GFAdd(vModSyndromes.data[iM],
                GFMpy(vSyndromes.data[iS], pErasures.data[iP]));
            iS++;}
        iS -= pErasures.size-1;}

genms0:
    printf("vModSyndromes:  ");
    for(iM = 0; iM < vGenRoots.size-vModSyndromes.size; iM++)
        printf("   ");
    ShowVector(&vModSyndromes);
    printf("pErasures:      ");
    ShowVector(&pErasures);
}

/*----------------------------------------------------------------------*/
/*      GenmModSyndromes        generate matrix for GenpErrorsM         */
/*                                                                      */
/*      In matrix, modified syndromes scroll down diagonally right      */
/*      to left. The inverted version of this matrix used to generate   */
/*      unknown polynomial.                                             */
/*----------------------------------------------------------------------*/
static void     GenmModSyndromes(void)
{
int     iM;                             /* VModSyn.. index */
int     iMM;                            /* MModSyn.. index */
int     iMMSize;
int     i;

    mModSyndromes.nrows = vErrLocators.size;
    mModSyndromes.ncols = mModSyndromes.nrows;
    iMMSize = mModSyndromes.ncols * mModSyndromes.nrows;

    iM = 0;
    for(iMM = 0; iMM < iMMSize;){
        for(i = 0; i < mModSyndromes.ncols; i++){
            mModSyndromes.data[iMM] = vModSyndromes.data[iM];
            iMM++;
            iM++;}
        iM -= mModSyndromes.ncols-1;}
}

/*----------------------------------------------------------------------*/
/*      GenpErrorsM     generate pErrors via matrix algorithm           */
/*                                                                      */
/*      Generate a matrix from vModSyndromes and invert it              */
/*      If inversion fails, reduce error count (until 0) and repeat     */
/*      pErrors = inverted matrix times higher order VModSyndromes      */
/*----------------------------------------------------------------------*/
static void GenpErrorsM(void)
{
int     iP;                             /* index to PErrors */
int     iM;                             /* index to VModSyndromes */
int     iI;                             /* index to IMMod... */
int     i;

    vErrLocators.size = vModSyndromes.size/2; /* init # errors */

    while(vErrLocators.size){           /* while # errors != 0 */
        GenmModSyndromes();             /* generate matrix */
#if DISPLAYM
        printf("mModSyndromes\n");
        ShowMatrix(&mModSyndromes);
#endif

        if(MatrixInv(&imModSyndromes,   /* invert matrix */
                     &mModSyndromes)){
            goto gnp1;}                 /* if fail try 1 less error */
#if DISPLAYM
        printf("imModSyndromes\n");
        ShowMatrix(&imModSyndromes);
#endif

        iM = imModSyndromes.nrows;      /* generate pErrors */
        pErrors.size = iM+1;
        pErrors.data[0] = 1;
        iI = 0;
        for(iP = imModSyndromes.nrows; iP; iP--){
            pErrors.data[iP] = 0;
            for(i = iI+imModSyndromes.ncols; iI < i;){
                pErrors.data[iP] = GFSub(pErrors.data[iP],
                    GFMpy(imModSyndromes.data[iI], vModSyndromes.data[iM]));
                iI++;
                iM++;}
            iM -= imModSyndromes.ncols;}

        if(Poly2Root(&vErrLocators, &pErrors)){ /* solve error poly */
gnp1:       vErrLocators.size--;                /* if fail try 1 less err */
            continue;}
        break;}

    if(!vErrLocators.size){             /* if no errors */
        pErrors.size = 1;               /*   handle no root case */
        pErrors.data[0] = 1;}

    printf("pErrors (M):    ");
    ShowVector(&pErrors);
}

/*----------------------------------------------------------------------*/
/*      GenpErrorsE     generate pErrors via Euclid division algorithm  */
/*----------------------------------------------------------------------*/
static void GenpErrorsE(void)
{
/* R[] is msb first | A[] is msb last (reversed) */
EUCLID *pED;                            /* R[i-2] | A[i-1] */
EUCLID *pER;                            /* R[i-1] | A[i-2] */
EUCLID *pET;                            /* temp */
int     i, j;
BYTE    bME;                            /* max errors possible */
BYTE    bQuot;                          /* quotient */

/*      E0 = initial ED: E0.R[-1] = x^MAXERR, E0.A[0] = 1 */
    E0.size = vModSyndromes.size+2;
    E0.indx = vModSyndromes.size+1;
    E0.data[0] = 1;
    memset(&E0.data[1], 0, vModSyndromes.size*sizeof(BYTE));
    E0.data[E0.indx] = 1;
    pED = &E0;

/*      E1 = initial ER: E1.R[0] = syndrome polynomial, E1.A[-1] = 0 */
    E1.size = vModSyndromes.size+2;
    E1.indx = vModSyndromes.size+1;
    E1.data[0] = 0;
    for(i = 1; i < E1.indx; i++){
        E1.data[i] = vModSyndromes.data[vModSyndromes.size-i];}
    E1.data[E1.indx] = 0;
    pER = &E1;

/*      init bME */
    bME = vModSyndromes.size/2;

/*      Euclid algorithm */

    while(1){                           /* while degree ER.R > max errors */ 
#if DISPLAYE
        printf("ED: ");
        ShowEuclid(pED);
        printf("ER: ");
        ShowEuclid(pER);
#endif
        while((pER->data[0] == 0) &&    /* shift dvsr left until msb!=0 */
              (pER->indx != 0)){        /*  or fully shifted left */
            pER->indx--;
            memcpy(&pER->data[0], &pER->data[1], (pER->size-1)*sizeof(BYTE));
            pER->data[pER->size-1] = 0;}

        if(pER->indx <= bME){           /* if degree ER.R[] <= bME, break */
            break;}

        while(1){                       /* while more sub-steps */
            if(pED->data[0]){           /*   if ED.R[] msb!=0, update ED, ER */
                bQuot = GFDiv(pED->data[0], pER->data[0]); /* Q=ED.R[msb]/ER.R[msb] */
                for(i = 0; i < pER->indx; i++){            /* ED.R[]=ED.R[]-Q*ER.R[] */
                    pED->data[i] = GFSub(pED->data[i], GFMpy(bQuot, pER->data[i]));}
                for(i = pED->indx; i < pER->size; i++){    /* ER.A[]=ER.A[]-Q*ED.A[] */
                    pER->data[i] = GFSub(pER->data[i], GFMpy(bQuot, pED->data[i]));}}
            if(pED->indx == pER->indx){ /*   if sub-steps done, break */
                break;}
            pED->indx--;                /*   shift ED left */
            memcpy(&pED->data[0], &pED->data[1], (pED->size-1)*sizeof(BYTE));
            pED->data[pED->size-1] = 0;}

        pET = pER;                      /* swap ED, ER */
        pER = pED;
        pED = pET;}

    pErrors.size = pED->size-pED->indx; /* set pErrors.size */

    if((pER->indx) >= pErrors.size){    /*  if degree ER.R too high */
        printf("GenpErrorsE remainder.size >= errors.size\n");
        goto fail0;}

    j = pErrors.size - 1;       /* right shift ER if Omega has leading zeroes */
    while(pER->indx < j){
        pER->indx++;
        for(i = pER->size-1; i;){
            i--;
            pER->data[i+1] = pER->data[i];}
        pER->data[0] = 0;}
#if DISPLAYE
    printf("EX: ");
    ShowEuclid(pER);
#endif

/*      pErrors = ED.A[] without unreversing = Lambda reversed */
    j = pED->indx;
    for(i = 0; i < pErrors.size; i++){
        pErrors.data[i] = pED->data[j];
        j++;}

#if DISPLAYE
    printf("pErrors (e):    ");
    ShowVector(&pErrors);
#endif

/*      Make most significant coef pErrors == 1 (divide by it) */
    bQuot = pErrors.data[0];
    if(bQuot == 0){
        printf("GenpErrorsE most sig coef of pErrors == 0\n");
        pLambda.size = 1;
        pLambda.data[0] = 1;
        goto fail0;}
    for(i = 0; i < pErrors.size; i++){
        pErrors.data[i] = GFDiv(pErrors.data[i], bQuot);}
    printf("pErrors (E):    ");
    ShowVector(&pErrors);

/*      if no erasures   */
/*       pOmega = ER.R[] */

    if(vErsLocators.size == 0){
        pOmega.size = pER->indx;
        for(i = 0; i < pOmega.size; i++){
            pOmega.data[i] = pER->data[i];}
        printf("pOmega  (e):    ");
        ShowVector(&pOmega);
        for(i = 0; i < pOmega.size; i++){
            pOmega.data[i] = GFDiv(pER->data[i], bQuot);}
        printf("pOmega  (E):    ");
        ShowVector(&pOmega);}

/*      check poly with all syndromes */

    bQuot = 0;
    j = pErrors.size;
    while(j <= vModSyndromes.size){
        for(i = 0; i < pErrors.size; i++){
            j--;
            bQuot = GFAdd(bQuot, GFMpy(pErrors.data[i], vModSyndromes.data[j]));}
        if(bQuot != 0){
            printf("GenpErrorsE syndrome check detected error\n");
            goto fail0;}
        j += pErrors.size+1;}

/*      Find roots of pErrors (if fail, then set for no roots) */

    if(Poly2Root(&vErrLocators, &pErrors)){ /* solve error poly */
        printf("GenpErrorsE poly2root failed \n");
fail0:
        pErrors.size = 1;                   /* handle no root case */
        pErrors.data[0] = 1;
        pOmega.size = 0;
        vErrLocators.size = 0;}
}

/*----------------------------------------------------------------------*/
/*      GenpErrorsB     generate pErrors via Berklekamp Massey          */
/*      note poly most signifcant index == 0                            */
/*----------------------------------------------------------------------*/
static void GenpErrorsB(void)
{
BYTE i, j, n;
BYTE L, m;
BYTE b, d;
BYTE db;

    b = 1;                          /* discrepancy when L last updated */
    L = 0;                          /* number of errors */
    m = 1;                          /* # iterations since L last updated */
    vB.size    = 1;
    vB.data[0] = 1;
    vC.size    = 1;
    vC.data[0] = 1;

    for(n = 0; n < vModSyndromes.size; n++){
        d = vModSyndromes.data[n];         /* calculate discrepancy */
        for(i = 1; i <= L; i++){
            d = GFAdd(d, GFMpy(vC.data[(vC.size - 1)- i], vModSyndromes.data[n-i]));}
        if(d == 0){                     /* if 0 increment m, continue */
            m += 1;
            continue;}
        vT.size = vC.size;              /* vT = vC */
        memcpy(vT.data, vC.data, vC.size);
        db = GFDiv(d,b);                /* db = (d/b) */
        vBx.size = vB.size+m;           /* Bx = x^m B */
        memcpy(vBx.data, vB.data, vB.size);
        memset(&vBx.data[vB.size], 0, m);
        for(i = 0; i < vBx.size; i++){  /* Bx *= db */
            vBx.data[i] = GFMpy(vBx.data[i], db);}
        j = vBx.size - vC.size;         /* right shift vBx or vC */
        if(((char)j) > 0){
            for(i = vBx.size; i > j; ){
                i--;
                vC.data[i] = vC.data[i-j];}
            memset(vC.data, 0, j);
            vC.size += j;}
        else if(((char)j) < 0){
            j = -j;
            for(i = vC.size; i > j; ){
                i--;
                vBx.data[i] = vBx.data[i-j];}
            memset(vBx.data, 0, j);
            vBx.size += j;}
        for(i = 0; i < vC.size; i++){   /* C -= Bx */
            vC.data[i] = GFSub(vC.data[i], vBx.data[i]);}
        if(n < 2*L){                    /* if L not increasing */
            m += 1;
            continue;}
        vB.size = vT.size;              /*   B = T */
        memcpy(vB.data, vT.data, vT.size);
        L = n + 1 - L;                  /* update L */
        b = d;                          /* update b */
        m = 1;}                         /* reset m */
}

/*----------------------------------------------------------------------*/
/*      MergeLocators           merge ers and err locators              */
/*----------------------------------------------------------------------*/
static void MergeLocators(void)
{
BYTE    bLcr;                           /* locator */
BYTE    bLoc;                           /* location */
int     xErs;                           /* indexes */
int     xErr;
int     xLoc;

    vLocators.size   = vErsLocators.size+vErrLocators.size;
    vLocations.size  = vLocators.size;
    vOffsets.size    = vLocators.size;

    bLcr = 1;                           /* set up */
    bLoc = 0;

    for(xErs = xErr = xLoc = 0; xLoc < vLocators.size;){

        if(xErs < vErsLocators.size &&
          bLcr == vErsLocators.data[xErs]){
            vLocators.data[xLoc] = bLcr;
            vLocations.data[xLoc]  = bLoc;
            xLoc++;
            xErs++;}
        if(xErr < vErrLocators.size &&
          bLcr == vErrLocators.data[xErr]){
            vLocators.data[xLoc] = bLcr;
            vLocations.data[xLoc]  = bLoc;
            xLoc++;
            xErr++;}

        bLcr = GFMpy(bLcr, bAlpha);      /* try next locator */
        bLoc++;}

    for(xLoc = 0; xLoc < vLocators.size; xLoc++){
        vOffsets.data[xLoc] = vData.size-1-vLocations.data[xLoc];}

    Root2Poly(&pLocators, &vLocators);  /* generate poly */
}

/*----------------------------------------------------------------------*/
/*      GenLambda                                                       */
/*----------------------------------------------------------------------*/
static void GenLambda(void)
{
BYTE i, j;
/*      pLambda = reverse of pLocators */

    pLambda.size = pLocators.size;
    for(i = 0; i < pLocators.size; i++){
        pLambda.data[i] = pLocators.data[(pLocators.size-1)-i];}

/*      generate vDLamda from pLambda  (copy odd terms)   */
/*          example: derivative of a x^3 + b x^2 + cx + d */
/*          (a+a+a) x^2 + (b+b) x + (c)                   */
/*                a x^2 +     0 x +  c                    */

    pDLambda.size = pLambda.size / 2;
    j = pDLambda.size-1;
    for(i = pLambda.size-2; i < (BYTE)0xfe;){
        pDLambda.data[j] = pLambda.data[i];
        j -= 1;
        i -= 2;}
}

/*----------------------------------------------------------------------*/
/*      GenOffsets  vInvLocators, vLocations, vOffsets                  */
/*----------------------------------------------------------------------*/
static void GenOffsets(void)
{
BYTE i;

    vInvLocators.size = vLocators.size;
    for(i = 0; i < vLocators.size; i++){
        vInvLocators.data[i] = GFDiv(1, vLocators.data[i]);}

    vLocations.size = vLocators.size;
    for(i = 0; i < vLocators.size; i++){
        vLocations.data[i] = abLog[vLocators.data[i]];}

    vOffsets.size = vLocators.size;
    for(i = 0; i < vLocators.size; i++){
        vOffsets.data[i] = vData.size-1-vLocations.data[i];}
}

/*----------------------------------------------------------------------*/
/*      GenmLocators    Generate Locator matrix (for error values)      */
/*                                                                      */
/*      Matrix is composed of vGenRoots powered to vLocations.          */
/*      The inverse of this matrix is used to produce error values.     */
/*      This method normally used for erasure only correction.          */
/*----------------------------------------------------------------------*/
static void GenmLocators(void)
{
BYTE    *p0;                            /* ptr to matrix */
BYTE    i, j;

    mLocators.nrows = vLocators.size;
    mLocators.ncols = vLocators.size;

    p0  = mLocators.data;               /* set up */

    for(j = 0; j < mLocators.nrows; j++){       /* do matrix */
        for(i = 0; i < mLocators.ncols; i++){
            *p0 = GFPow(vGenRoots.data[j],  vLocations.data[i]);
            p0++;}}
}

/*----------------------------------------------------------------------*/
/*      GenErrValues    generate error values                           */
/*                                                                      */
/*      Multiplying inverse locator matrix by the first few syndromes   */
/*      produces error values.                                          */
/*----------------------------------------------------------------------*/
static void GenErrValues(void)
{
BYTE    bE;                             /* index to vErrValues */
BYTE    bS;                             /* index to vSyndromes */
BYTE    bI;                             /* index to imLocators */
BYTE    i;

    bI = 0;
    vErrValues.size = imLocators.nrows;

    for(bE = 0; bE < vErrValues.size; bE++){
        vErrValues.data[bE] = 0;
        bS = 0;
        for(i = 0; i < imLocators.ncols; i++){
            vErrValues.data[bE] = GFAdd(vErrValues.data[bE],
                GFMpy(imLocators.data[bI], vSyndromes.data[bS]));
            bI++;
            bS++;}}
}

/*----------------------------------------------------------------------*/
/*      GenAltInvLoc    Generate alternate inverted locator matrix      */
/*                      and alternate error values                      */
/*                                                                      */
/*      Uses pLocators and vLocators to create imLocators equivalent.   */
/*      This method normally used for erasure only correction.          */
/*----------------------------------------------------------------------*/
static void GenAltInvLoc(void)
{
BYTE    *pm;                            /* ptr to matrix */
BYTE    b0, b1;                         /* temps */
BYTE    i, j;

    mAltInvLoc.nrows = vLocators.size;
    mAltInvLoc.ncols = vLocators.size;
    vAltErrValues.size = vLocators.size;

    pm = mAltInvLoc.data;               /* set up */
    for(j = 0; j < mAltInvLoc.nrows; j++){

/*      Add a row: reversed poly divide (pLocators)/(1X - vLocators[j])   */
/*      This is equivalent to reversed order coefs of:                  */
/*      Root2Poly(row, vLocators with "jth" element dropped).           */

        pm += vLocators.size;
        b0 = 0;
        for(i = 0; i < mAltInvLoc.ncols; i++){
            *--pm = b0 = GFAdd(pLocators.data[i], GFMpy(b0, vLocators.data[j]));}

/*      b0 = sum of row[i]*Locator[j]**i = divisor for next step */

        b0 = 0;
        for(i = 0; i < mAltInvLoc.ncols; i++){
            b0 = GFAdd(b0, GFMpy(*pm, GFPow(vLocators.data[j],i)));
            pm++;}
        pm -= vLocators.size;

/*      adjust b0 for first root of vGenRoots */

        b0 = GFMpy(b0,GFPow(vGenRoots.data[0], vLocations.data[j]));

/*      divide the row by b0 */
/*      error value bS summed up in b1. The 3 lines below       */
/*      using b1 are not needed for matrix generation.          */

        b1 = 0;
        for(i = 0; i < mAltInvLoc.ncols; i++){
            *pm = GFDiv(*pm, b0);
            b1 = GFAdd(b1, GFMpy(*pm, vSyndromes.data[i]));
            pm++;}
        vAltErrValues.data[j] = b1;
    }
}

/*----------------------------------------------------------------------*/
/*      GenOmega        generate Omega                                  */
/*----------------------------------------------------------------------*/
static void GenOmega(void)
{
int iO;                                 /* index to pOmega */
int iL;                                 /* index to pLambda */
int iS;                                 /* index to vSyndromes */
    pOmega.size = vLocators.size;
    for(iO = vLocators.size-1; iO >= 0; iO--){
        iS = vLocators.size-1-iO;
        pOmega.data[iO] = vSyndromes.data[iS--];
        for(iL = vLocators.size-1; iL > iO; iS, iL--){
            pOmega.data[iO] = GFAdd(pOmega.data[iO],
                    GFMpy(vSyndromes.data[iS--], pLambda.data[iL]));}}
}

/*----------------------------------------------------------------------*/
/*      GenForneyErr    generate error values using Forney              */
/*----------------------------------------------------------------------*/
static void GenForneyErr(void)
{
BYTE    i,j;
BYTE    bDvnd;
BYTE    bDvsr;

    vForney.size = pOmega.size;
    memset(vForney.data, 0, vForney.size);
    for(j = 0; j < vForney.size; j++){
        bDvsr = bDvnd = 0;
        for(i = pOmega.size; i;){
            i--;
            bDvnd = GFAdd(bDvnd, GFMpy(pOmega.data[i],
                    GFPow(vInvLocators.data[j], (BYTE)((pOmega.size-1)-i))));}
        bDvnd = GFMpy(bDvnd, GFPow(vLocators.data[j], bForneyC));
        for(i = pDLambda.size; i;){
            i--;
            bDvsr = GFAdd(bDvsr, GFMpy(pDLambda.data[i],
                    GFPow(vInvLocators.data[j], GFMpy(2, (BYTE)((pDLambda.size-1)-i)))));}
        if(bDvsr == 0){
            printf("GenForneyErr divide by zero\n");
            vForney.size = 0;
            return;}
        vForney.data[j] = GFSub(0, GFDiv(bDvnd, bDvsr));}
}

/*----------------------------------------------------------------------*/
/*      FixvData        Fix up vData                                    */
/*----------------------------------------------------------------------*/
static void FixvData(void)
{
BYTE    i;

    for(i = 0; i < vOffsets.size; i++){
        vData.data[vOffsets.data[i]] = GFSub(
            vData.data[vOffsets.data[i]],
            vErrValues.data[i]);}
}

/*----------------------------------------------------------------------*/
/*      Poly2Root(pVDst, pPSrc)         find roots of poly              */
/*----------------------------------------------------------------------*/
static int Poly2Root(VECTOR *pVDst, VECTOR *pPSrc)
{
BYTE    bLcr;                           /* current locator */
BYTE    bSum;                           /* current sum */
BYTE    bV;                             /* index to pVDst */
BYTE    i,j;

    pVDst->size = pPSrc->size-1;        /* set dest size */

    if(!pVDst->size)                    /* exit if null */
        return(0);

    bV   = 0;
    bLcr = 1;
    for(j = 0; j < vData.size;  j++){
        bSum = 0;                       /* sum up terms */
        for(i = 0; i < pPSrc->size; i++){
            bSum = GFMpy(bSum, bLcr);
            bSum = GFAdd(bSum, pPSrc->data[i]);}

        if(!bSum){                      /* if a root */
            if(bV == pVDst->size){      /*    exit if too many roots */
                return(1);}
            pVDst->data[bV] = bLcr;     /*    add locator */
            bV++;}

        bLcr = GFMpy(bLcr, bAlpha);}    /* set up next higher Alpha */

    if(bV != pVDst->size)               /* exit if not enough roots */
        return(1);

    return(0);
}

/*----------------------------------------------------------------------*/
/*      Root2Poly(pPDst, pVSrc)         convert roots into polynomial   */
/*----------------------------------------------------------------------*/
static void Root2Poly(VECTOR *pPDst, VECTOR *pVSrc)
{
BYTE i, j;

    pPDst->size = pVSrc->size+1;
    pPDst->data[0] = 1;
    memset(&pPDst->data[1], 0, pVSrc->size);
    for(j = 0; j < pVSrc->size; j++){
        for(i = j; i < (BYTE)0xff; i--){
            pPDst->data[i+1] = GFSub(pPDst->data[i+1],
                    GFMpy(pPDst->data[i], pVSrc->data[j]));}}
}

/*----------------------------------------------------------------------*/
/*      MatrixInv(pMDst, pMSrc) invert matrix                           */
/*      assumes square matrix                                           */
/*----------------------------------------------------------------------*/
static int MatrixInv(MATRIX *pMDst, MATRIX *pMSrc)
{
BYTE    *p0;                            /* generic pointers */
BYTE    *p1;
BYTE    *p2;
BYTE    bMod;                           /* row modifier */
BYTE    bTemp;
BYTE    bNCol;                          /* # columns */
BYTE    bNCol2;                         /* 2 * # columns */
BYTE    bNCol21;                        /* 1+2*# columns */
BYTE    i, j;

    memset(abId, 0, sizeof(abId));      /* set up abId */
    abId[MAXPAR] =  1;

    bNCol   = pMSrc->ncols;             /* set size related stuff */
    bNCol2  = bNCol+bNCol;
    bNCol21 = bNCol2+1;

    pMDst->nrows = bNCol;               /* set destination size */
    pMDst->ncols = bNCol2;              /*   for display */

/*      generate double-width augmented matrix */
/*      left side is copy of pMSrc, right side is Identity Matrix */

    p0  = pMSrc->data;
    p1  = pMDst->data;
    for(i = 0; i < bNCol; i++){
        memcpy(p1,  p0, bNCol); /* copy a row of data */
        p0  += bNCol;
        p1  += bNCol;
        memcpy(p1,  &abId[MAXPAR-i], bNCol); /* add ID  part */
        p1  += bNCol;}

/*      normalize according to left side */
/*      results in inverse matrix in right size */

    #if DISPLAYI
        printf("start\n");
        ShowMatrix(pMDst);
    #endif

    for(j = 0; j < bNCol; j++){         /* working at [j][j] */

/*      find 1st non-zero in current column */

        p0  = pMDst->data+j*bNCol21;    /* p0 = starting ptr */
        p1  = p0;                       /* p1 = scanning ptr */
        p2  = pMDst->data+bNCol2*bNCol; /* p2 = end of  matrix */
        while(0 == *p1){
            p1  += bNCol2;
            if(p1 >= p2){               /* return if bad matrix */
                return(1);}}

        bMod = *p1;                     /* set divisor for row */

/*      swap rows if needed */

        if(p0 != p1){
            p0  -= j;                   /* backup to beg of rows */
            p1  -= j;
            for(p2  = p0+bNCol2; p0 != p2; p0++, p1++){
                bTemp = *p0;
                *p0 = *p1;
                *p1 = bTemp;}
            #if DISPLAYI
                printf("swapped rows\n");
                ShowMatrix(pMDst);
            #endif
            ;}


/*      divide row to produce a one */

        p0  = pMDst->data+j*bNCol2;     /* p0 = ptr to  start of row */
        for(p2  = p0+bNCol2; p0 != p2; p0++){
            *p0 = GFDiv(*p0, bMod);}

        #if DISPLAYI
            printf("divided row\n");
            ShowMatrix(pMDst);
        #endif

/*      subtract multiple of this row from other rows */
/*      to create a column of zeroes */
/*      (except for this row,column) */

        for(i = 0; i < bNCol; i++){
            if(i == j)                  /* skip if current row */
                continue;
            p0  = pMDst->data+j*bNCol2; /* p0 = ptr to  current row */
            p1  = pMDst->data+i*bNCol2; /* p1 = ptr to  target row */
            bMod = *(p1+j);             /* bMod = current row mpyr */
            for(p2  = p0+bNCol2; p0 != p2; p0++, p1++)
                *p1 = GFSub(*p1, GFMpy(*p0, bMod));}
#if DISPLAYI
        printf("zeroed columns\n");
        ShowMatrix(pMDst);
#endif
        ;}

/*      now copy right side of matrix to left side */

    p0  = pMDst->data;                  /* p0 = Dst ptr */
    p1  = p0+bNCol;
    for(j = 0; j < bNCol; j++){
        p2  = p0+bNCol;
        while(p0 != p2){
            *p0++ = *p1++;}
        p1  += bNCol;}

    pMDst->ncols = bNCol;               /* set proper ncols */
    return(0);
}

/*----------------------------------------------------------------------*/
/*      GFAdd(b0, b1)           b0+b1                                   */
/*----------------------------------------------------------------------*/
static BYTE GFAdd(BYTE b0, BYTE b1)
{
    return(b0^b1);
}

/*----------------------------------------------------------------------*/
/*      GFSub(b0, b1)           b0-b1                                   */
/*----------------------------------------------------------------------*/
static BYTE GFSub(BYTE b0, BYTE b1)
{
    return(b0^b1);
}
/*----------------------------------------------------------------------*/
/*      GFMpy(b0, b1)           b0*b1       using logs                  */
/*----------------------------------------------------------------------*/
static BYTE GFMpy(BYTE byt0, BYTE byt1)
{
    if(byt0 == 0 || byt1 == 0)
        return(0);

    return(abExp[(int)abLog[byt0]+(int)abLog[byt1]]);
}

/*----------------------------------------------------------------------*/
/*      GFDiv(b0, b1)           b0/b1                                   */
/*----------------------------------------------------------------------*/
static BYTE GFDiv(BYTE b0, BYTE b1)
{
    if(b1 == 0){
        printf("divide by zero\n");
        return(0);}
    if(b0 == 0)
        return(0);
    return(abExp[255+(int)abLog[b0]-(int)abLog[b1]]);
}

/*----------------------------------------------------------------------*/
/*      GFPow(b0, b1)           b0^b1                                   */
/*----------------------------------------------------------------------*/
static BYTE GFPow(BYTE b0, BYTE b1)
{
BYTE b;
    b = 1;
    while(b1){
        if(b1&1)
            b = GFMpy(b, b0);
        b0 = GFMpy(b0, b0);
        b1 >>= 1;}
    return(b);
}

/*----------------------------------------------------------------------*/
/*      GFMpy0(b0,b1)           b0*b1       using low level math        */
/*----------------------------------------------------------------------*/
static BYTE GFMpy0(BYTE b0, BYTE b1)
{
int i;
int product;
    product = 0;
    for(i = 0; i < 8; i++){
        product <<= 1;
        if(product & 0x100){
            product ^= iGF;}
        if(b0 & 0x80u){
            product ^= b1;}
        b0 <<= 1;}
    return((BYTE)product);
}

/*----------------------------------------------------------------------*/
/*      InitGF  Initialize Galios Stuff                                 */
/*----------------------------------------------------------------------*/
static void InitGF(void)
{
BYTE b;
int i;
int j;

    b = 1;
    for(i = 0; i < 512; i++){           /* init abExp[] */
        abExp[i] = b;
        b = GFMpy0(b, bAlpha);}

    abLog[0] = 0xff;                    /* init abLog[] */
    for(i = 0; i < 255; i++){
        abLog[abExp[i]] = i;}

/*  b = first consective root (fcr)                   */
    if(!fSelfRcp){                      /* if not self-reciprocal poly */
        b = bFCR;                       /* normal narrow sense fcr = Alpha */
    } else {                            /* else self reciprocal poly */
        if(vGenRoots.size&1){           /*   set starting root       */
            j = 255-(vGenRoots.size>>1);}
        else{
            j = 128-(vGenRoots.size>>1);}
        b = GFPow(bAlpha, (BYTE)j);
        bFCR = b;}

    for(i = 0; i < vGenRoots.size; i++){
        vGenRoots.data[i] = b;
        b = GFMpy(b, bAlpha);}

    Root2Poly(&pGenPoly, &vGenRoots);   /* convert to poly */

    j = bLogFCR = abLog[vGenRoots.data[0]];
    j = 1 - j;                          /* set forney correction factor */
    if(j < 0)j += 255;
    bForneyC = (BYTE)j;
}

/*----------------------------------------------------------------------*/
/*      ShowEuclid                                                      */
/*----------------------------------------------------------------------*/
static void ShowEuclid(EUCLID *pESrc)
{
BYTE    i;
    for(i = 0; i < pESrc->indx; i++){
        printf(" %02x", pESrc->data[i]);}
    printf("|");
    for( ; i < pESrc->size; i++){
        printf("%02x ", pESrc->data[i]);}
    printf("\n");
}

/*----------------------------------------------------------------------*/
/*      ShowVector                                                      */
/*----------------------------------------------------------------------*/
static void ShowVector(VECTOR *pVSrc)
{
BYTE    i;
    for(i = 0; i < pVSrc->size; ){
        printf(" %02x", pVSrc->data[i]);
        i++;
        if(0 == (i&0xf)){
            printf("\n");}}
    printf("\n");
}

/*----------------------------------------------------------------------*/
/*      ShowMatrix                                                      */
/*----------------------------------------------------------------------*/
static void ShowMatrix(MATRIX *pMSrc)
{
BYTE    bM;
BYTE    i, j;
    bM = 0;
    for(i = 0; i < pMSrc->nrows; i++){
        for(j = 0; j < pMSrc->ncols; j++){
            printf(" %02x", pMSrc->data[bM++]);}
        printf("\n");}
}

/*----------------------------------------------------------------------*/
/*      Conrs get console response - emulates _cgets()                  */
/*----------------------------------------------------------------------*/
static int Conrs(void)
{
int i;
    memset(abUser, 0, sizeof(abUser));  /* get a line */
    abUser[0] = sizeof(abUser)-2;
    fgets(abUser+2, sizeof(abUser)-2, stdin);
    abUser[1] = (BYTE)(strlen(&abUser[2])-1);
    i = abUser[1];
    abUser[2+i] = 0;
    return(i);
}

/*----------------------------------------------------------------------*/
/*      main                                                            */
/*----------------------------------------------------------------------*/
int main()
{
    DoUser();
    return(0);
}

/*----------------------------------------------------------------------*/
/*      DoUser  do user stuff                                           */
/*----------------------------------------------------------------------*/

static void DoUser(void)
{
int i, j, k;
int iU0, iU1;

    printf("eccdemo8 1.7\n");
/*      select Galios Field */
    while(1){
        printf("Galios Field Generators and Alpha's\n");
        printf("   ");
        for(i = 0; i < 10; i++){
            printf("%7d", i);}
        printf("\n");
        iU0 = 0;
        for(j = 0; j < 3; j++){
            printf("%02d:", iU0/2);
            for(i = 0; i < 10; i++){
                printf(" %03x %02x", aiGA[iU0], aiGA[iU0+1]);
                iU0 += 2;}
            printf("\n");}
        printf("Select Galios Field Generator (0-29):\n");
        if(!Conrs())return;
        sscanf(&abUser[2], "%d", &iU0);
        if(iU0 > 29)continue;
        iU0  <<= 1;
        iGF    = aiGA[iU0];
        bAlpha = (BYTE)(aiGA[iU0+1]);
        printf(    "    %03x %02x\n", iGF, bAlpha);
        break;}

/*      Select normal or self-reciprocal poly */
    {
        fSelfRcp = 0;
        printf("Use self-reciprocal poly? (Y/N):\n");
        if(!Conrs())return;
        if('Y' == (0x5f&abUser[2])){
            fSelfRcp = 1;}
    }

/*      Select first consecutive root */
    if(fSelfRcp == 0){
        printf("Enter first consecutive root == 1 or %1x:\n", bAlpha);
        if(!Conrs())return;
        if('1' == abUser[2]){
            bFCR = 1;}
        else{
            bFCR = bAlpha;
        }
    }

    while(1){
        printf("Select # of parity bytes (1 -> %d):\n", MAXPAR);
        if(!Conrs())return;
        sscanf(&abUser[2], "%d", &iU0);
        if(iU0 == 0 || iU0 > MAXPAR)continue;
        vGenRoots.size = iU0;
        iU1 = 255-vGenRoots.size;
        InitGF();
        break;}

    while(1){
        printf("Select # of data bytes (1 -> %d):\n", iU1);
        if(!Conrs())return;
        sscanf(&abUser[2], "%d", &iU0);
        if(iU0 == 0 || iU0 > iU1)continue;
        iNumData = iU0;
        vData.size = iU0+vGenRoots.size;
        break;}

DoUser0:

    printf(    "vGenRoots:      ");
    ShowVector(&vGenRoots);
    printf(    "pGenPoly:       ");
    ShowVector(&pGenPoly);
    memset(vData.data, 0, vData.size);  /* reset vData */

DoUser1:

    vErsf.size = vData.size;            /* reset vErsf */
    memset(vErsf.data, 0, vErsf.size*sizeof(BYTE));

DoUser2:

/*      display data */

    printf("vErsf:\n");
    ShowVector(&vErsf);
    printf("vData:\n");
    ShowVector(&vData);

DoUser3:

/*      get user command */
    printf("Enter Change Ptr Zero Ncode Fix Quit: ");
    if(0 == Conrs())
        goto DoUser3;
    printf("\n");
    i = abUser[2] & 0x5f;
    switch(i){
      case 'E':                         /* enter data */
        for(i = 0; i < vData.size; i++){
            printf("%2d = ", i);
            if(!Conrs())break;
            sscanf(&abUser[2], "%x", &iU0);
            vData.data[i] = iU0;}
        break;
      case 'C':                         /* change data */
        k = 0;
DoUserc:
        while(1){
            printf("Enter offset (0-%d): ", vData.size-1);
            if(!Conrs())break;
            sscanf(&abUser[2], "%d", &iU0);
            if(iU0 >= vData.size)continue;
            printf("Enter value: (0-ff): ");
            if(!Conrs())break;
            sscanf(&abUser[2], "%x", &iU1);
            vData.data[iU0] = iU1;
            vErsf.data[iU0] = k;}
        break;
      case 'P':                         /* pointers */
        k = 1;
        goto DoUserc;
      case 'Z':                         /* zero data array */
        goto DoUser0;
      case 'N':                         /* encode ecc bytes */
        Encode();
        memcpy(&vData.data[iNumData], vParities.data, vParities.size);
        break;
      case 'F':                         /* fix data */
        Decode();
        goto DoUser1;
        break;
      case 'Q':                         /* quit */
        return;
      default:
        break;}
    goto DoUser2;
    return;
}
