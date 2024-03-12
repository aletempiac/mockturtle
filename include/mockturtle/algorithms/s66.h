/*************************************************************
   Decomposition of 11-input functions into two 6-input boxes
**************************************************************/

// This stand-alone code is adapted from ABC (file "src/map/if/ifDec16.c")
// https://github.com/berkeley-abc/abc/blob/master/src/map/if/ifDec16.c

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <time.h>


#ifdef _WIN32
typedef unsigned __int64 word;   // 32-bit windows
#else
typedef long long unsigned word; // other platforms
#endif

#ifdef LIN
#define ABC_CONST(number) number ## ULL 
#else // LIN64 and windows
#define ABC_CONST(number) number
#endif


////////////////////////////////////////////////////////////////////////
///                        DECLARATIONS                              ///
////////////////////////////////////////////////////////////////////////

#define CLU_VAR_MAX  11
#define CLU_WRD_MAX  (1 << ((CLU_VAR_MAX)-6))
#define CLU_MEM_MAX  1000  // 1 GB
#define CLU_UNUSED   0xff

// decomposition
typedef struct If_Grp_t_ If_Grp_t;
struct If_Grp_t_
{
    char       nVars;
    char       nMyu;
    char       pVars[CLU_VAR_MAX];
};

// variable swapping code
static word PMasks[5][3] = {
    { ABC_CONST(0x9999999999999999), ABC_CONST(0x2222222222222222), ABC_CONST(0x4444444444444444) },
    { ABC_CONST(0xC3C3C3C3C3C3C3C3), ABC_CONST(0x0C0C0C0C0C0C0C0C), ABC_CONST(0x3030303030303030) },
    { ABC_CONST(0xF00FF00FF00FF00F), ABC_CONST(0x00F000F000F000F0), ABC_CONST(0x0F000F000F000F00) },
    { ABC_CONST(0xFF0000FFFF0000FF), ABC_CONST(0x0000FF000000FF00), ABC_CONST(0x00FF000000FF0000) },
    { ABC_CONST(0xFFFF00000000FFFF), ABC_CONST(0x00000000FFFF0000), ABC_CONST(0x0000FFFF00000000) }
};
// elementary truth tables
static word Truth6[6] = {
    ABC_CONST(0xAAAAAAAAAAAAAAAA),
    ABC_CONST(0xCCCCCCCCCCCCCCCC),
    ABC_CONST(0xF0F0F0F0F0F0F0F0),
    ABC_CONST(0xFF00FF00FF00FF00),
    ABC_CONST(0xFFFF0000FFFF0000),
    ABC_CONST(0xFFFFFFFF00000000)
};
static word Truths6Neg[6] = {
    ABC_CONST(0x5555555555555555),
    ABC_CONST(0x3333333333333333),
    ABC_CONST(0x0F0F0F0F0F0F0F0F),
    ABC_CONST(0x00FF00FF00FF00FF),
    ABC_CONST(0x0000FFFF0000FFFF),
    ABC_CONST(0x00000000FFFFFFFF)
};

static word TruthAll[CLU_VAR_MAX][CLU_WRD_MAX] = {{0}};


extern void Kit_DsdPrintFromTruth( unsigned * pTruth, int nVars );

//void Extra_PrintBinary( FILE * pFile, unsigned Sign[], int nBits ) {}

extern int If_CluSupportSize( word * t, int nVars );

int s_Count2 = 0;
int s_Count3 = 0;

static inline int      Abc_AbsInt( int a        )             { return a < 0 ? -a : a; }
static inline int      Abc_MaxInt( int a, int b )             { return a > b ?  a : b; }
static inline int      Abc_MinInt( int a, int b )             { return a < b ?  a : b; }

static inline int  Abc_TtWordNum( int nVars )     { return nVars <= 6 ? 1 : 1 << (nVars-6); }
static inline int  Abc_TtHexDigitNum( int nVars ) { return nVars <= 2 ? 1 : 1 << (nVars-2); }



static inline int Abc_TtHasVar( word * t, int nVars, int iVar )
{
    assert( iVar < nVars );
    assert( nVars > 6 );
    //if ( nVars <= 6 )
    //    return Abc_Tt6HasVar( t[0], iVar );
    if ( iVar < 6 )
    {
        int i, Shift = (1 << iVar);
        int nWords = Abc_TtWordNum( nVars );
        for ( i = 0; i < nWords; i++ )
            if ( ((t[i] >> Shift) & Truths6Neg[iVar]) != (t[i] & Truths6Neg[iVar]) )
                return 1;
        return 0;
    }
    else
    {
        int i, Step = (1 << (iVar - 6));
        word * tLimit = t + Abc_TtWordNum( nVars );
        for ( ; t < tLimit; t += 2*Step )
            for ( i = 0; i < Step; i++ )
                if ( t[i] != t[Step+i] )
                    return 1;
        return 0;
    }
}
static inline word Abc_Tt6Stretch( word t, int nVars )
{
    assert( nVars >= 0 );
    if ( nVars == 0 )
        nVars++, t = (t & 0x1) | ((t & 0x1) << 1);
    if ( nVars == 1 )
        nVars++, t = (t & 0x3) | ((t & 0x3) << 2);
    if ( nVars == 2 )
        nVars++, t = (t & 0xF) | ((t & 0xF) << 4);
    if ( nVars == 3 )
        nVars++, t = (t & 0xFF) | ((t & 0xFF) << 8);
    if ( nVars == 4 )
        nVars++, t = (t & 0xFFFF) | ((t & 0xFFFF) << 16);
    if ( nVars == 5 )
        nVars++, t = (t & 0xFFFFFFFF) | ((t & 0xFFFFFFFF) << 32);
    assert( nVars == 6 );
    return t;
}

////////////////////////////////////////////////////////////////////////
///                     FUNCTION DEFINITIONS                         ///
////////////////////////////////////////////////////////////////////////

static inline unsigned If_CluGrp2Uns( If_Grp_t * pG )
{
    char * pChar = (char *)pG;
    unsigned Res = 0;
    int i;
    for ( i = 0; i < 8; i++ )
        Res |= ((pChar[i] & 15) << (i << 2));
    return Res;
}

static inline void If_CluUns2Grp( unsigned Group, If_Grp_t * pG )
{
    char * pChar = (char *)pG;
    int i;
    for ( i = 0; i < 8; i++ )
        pChar[i] = ((Group >> (i << 2)) & 15);
}

unsigned int If_CluPrimeCudd( unsigned int p )
{
    int i,pn;

    p--;
    do {
        p++;
        if (p&1) {
        pn = 1;
        i = 3;
        while ((unsigned) (i * i) <= p) {
        if (p % i == 0) {
            pn = 0;
            break;
        }
        i += 2;
        }
    } else {
        pn = 0;
    }
    } while (!pn);
    return(p);

} /* end of Cudd_Prime */

// hash table
static inline int If_CluWordNum( int nVars )
{
    return nVars <= 6 ? 1 : 1 << (nVars-6);
}
static inline int If_CluCountOnes( word t )
{
    t =    (t & ABC_CONST(0x5555555555555555)) + ((t>> 1) & ABC_CONST(0x5555555555555555));
    t =    (t & ABC_CONST(0x3333333333333333)) + ((t>> 2) & ABC_CONST(0x3333333333333333));
    t =    (t & ABC_CONST(0x0F0F0F0F0F0F0F0F)) + ((t>> 4) & ABC_CONST(0x0F0F0F0F0F0F0F0F));
    t =    (t & ABC_CONST(0x00FF00FF00FF00FF)) + ((t>> 8) & ABC_CONST(0x00FF00FF00FF00FF));
    t =    (t & ABC_CONST(0x0000FFFF0000FFFF)) + ((t>>16) & ABC_CONST(0x0000FFFF0000FFFF));
    return (t & ABC_CONST(0x00000000FFFFFFFF)) +  (t>>32);
}

// variable permutation for large functions
static inline void If_CluClear( word * pIn, int nVars )
{
    int w, nWords = If_CluWordNum( nVars );
    for ( w = 0; w < nWords; w++ )
        pIn[w] = 0;
}
static inline void If_CluFill( word * pIn, int nVars )
{
    int w, nWords = If_CluWordNum( nVars );
    for ( w = 0; w < nWords; w++ )
        pIn[w] = ~(word)0;
}
static inline void If_CluCopy( word * pOut, word * pIn, int nVars )
{
    int w, nWords = If_CluWordNum( nVars );
    for ( w = 0; w < nWords; w++ )
        pOut[w] = pIn[w];
}
static inline int If_CluEqual( word * pOut, word * pIn, int nVars )
{
    int w, nWords = If_CluWordNum( nVars );
    for ( w = 0; w < nWords; w++ )
        if ( pOut[w] != pIn[w] )
            return 0;
    return 1;
}
static inline void If_CluAnd( word * pRes, word * pIn1, word * pIn2, int nVars )
{
    int w, nWords = If_CluWordNum( nVars );
    for ( w = 0; w < nWords; w++ )
        pRes[w] = pIn1[w] & pIn2[w];
}
static inline void If_CluSharp( word * pRes, word * pIn1, word * pIn2, int nVars )
{
    int w, nWords = If_CluWordNum( nVars );
    for ( w = 0; w < nWords; w++ )
        pRes[w] = pIn1[w] & ~pIn2[w];
} 
static inline void If_CluOr( word * pRes, word * pIn1, word * pIn2, int nVars )
{
    int w, nWords = If_CluWordNum( nVars );
    for ( w = 0; w < nWords; w++ )
        pRes[w] = pIn1[w] | pIn2[w];
}
static inline word If_CluAdjust( word t, int nVars )
{
    assert( nVars >= 0 && nVars <= 6 );
    if ( nVars == 6 )
        return t;
    t &= (((word)1) << (1 << nVars)) - 1;
    if ( nVars == 0 )
        t |= t << (1<<nVars++);
    if ( nVars == 1 )
        t |= t << (1<<nVars++);
    if ( nVars == 2 )
        t |= t << (1<<nVars++);
    if ( nVars == 3 )
        t |= t << (1<<nVars++);
    if ( nVars == 4 )
        t |= t << (1<<nVars++);
    if ( nVars == 5 )
        t |= t << (1<<nVars++);
    return t;
}
static inline void If_CluAdjustBig( word * pF, int nVarsCur, int nVarsMax )
{
    int v, nWords;
    if ( nVarsCur == nVarsMax )
        return;
    assert( nVarsCur < nVarsMax );
    for ( v = Abc_MaxInt( nVarsCur, 6 ); v < nVarsMax; v++ )
    {
        nWords = If_CluWordNum( v );
        If_CluCopy( pF + nWords, pF, v );
    }
}
static inline void If_CluSwapAdjacent( word * pOut, word * pIn, int iVar, int nVars )
{
    int i, k, nWords = If_CluWordNum( nVars );
    assert( iVar < nVars - 1 );
    if ( iVar < 5 )
    {
        int Shift = (1 << iVar);
        for ( i = 0; i < nWords; i++ )
            pOut[i] = (pIn[i] & PMasks[iVar][0]) | ((pIn[i] & PMasks[iVar][1]) << Shift) | ((pIn[i] & PMasks[iVar][2]) >> Shift);
    }
    else if ( iVar > 5 )
    {
        int Step = (1 << (iVar - 6));
        for ( k = 0; k < nWords; k += 4*Step )
        {
            for ( i = 0; i < Step; i++ )
                pOut[i] = pIn[i];
            for ( i = 0; i < Step; i++ )
                pOut[Step+i] = pIn[2*Step+i];
            for ( i = 0; i < Step; i++ )
                pOut[2*Step+i] = pIn[Step+i];
            for ( i = 0; i < Step; i++ )
                pOut[3*Step+i] = pIn[3*Step+i];
            pIn  += 4*Step;
            pOut += 4*Step;
        }
    }
    else // if ( iVar == 5 )
    {
        for ( i = 0; i < nWords; i += 2 )
        {
            pOut[i]   = (pIn[i]   & ABC_CONST(0x00000000FFFFFFFF)) | ((pIn[i+1] & ABC_CONST(0x00000000FFFFFFFF)) << 32);
            pOut[i+1] = (pIn[i+1] & ABC_CONST(0xFFFFFFFF00000000)) | ((pIn[i]   & ABC_CONST(0xFFFFFFFF00000000)) >> 32);
        }
    }
}

void If_CluChangePhase( word * pF, int nVars, int iVar )
{
    int nWords = If_CluWordNum( nVars );
    assert( iVar < nVars );
    if ( iVar < 6 )
    {
        int i, Shift = (1 << iVar);
        for ( i = 0; i < nWords; i++ )
            pF[i] = ((pF[i] & ~Truth6[iVar]) << Shift) | ((pF[i] & Truth6[iVar]) >> Shift);
    }
    else
    {
        word Temp;
        int i, k, Step = (1 << (iVar - 6));
        for ( k = 0; k < nWords; k += 2*Step )
        {
            for ( i = 0; i < Step; i++ )
            {
                Temp = pF[i];
                pF[i] = pF[Step+i];
                pF[Step+i] = Temp;
            }
            pF += 2*Step;
        }
    }
}
void If_CluCountOnesInCofs( word * pTruth, int nVars, int * pStore )
{
    int nWords = If_CluWordNum( nVars );
    int i, k, nOnes = 0, Limit = Abc_MinInt( nVars, 6 );
    memset( pStore, 0, sizeof(int) * 2 * nVars );
    // compute positive cofactors
    for ( k = 0; k < nWords; k++ )
        for ( i = 0; i < Limit; i++ )
            pStore[2*i+1] += If_CluCountOnes( pTruth[k] & Truth6[i] );
    if ( nVars > 6 )
    for ( k = 0; k < nWords; k++ )
        for ( i = 6; i < nVars; i++ )
            if ( k & (1 << (i-6)) )
                pStore[2*i+1] += If_CluCountOnes( pTruth[k] );
    // compute negative cofactors
    for ( k = 0; k < nWords; k++ )
        nOnes += If_CluCountOnes( pTruth[k] );
    for ( i = 0; i < nVars; i++ )
        pStore[2*i] = nOnes - pStore[2*i+1];
}
unsigned If_CluSemiCanonicize( word * pTruth, int nVars, int * pCanonPerm )
{
    word pFunc[CLU_WRD_MAX], * pIn = pTruth, * pOut = pFunc, * pTemp;
    int pStore[CLU_VAR_MAX*2];
    unsigned uCanonPhase = 0;
    int i, Temp, fChange, Counter = 0;
//Kit_DsdPrintFromTruth( (unsigned*)pTruth, nVars ); printf( "\n" );

    // collect signatures 
    If_CluCountOnesInCofs( pTruth, nVars, pStore );
    // canonicize phase
    for ( i = 0; i < nVars; i++ )
    {
        if ( pStore[2*i+0] <= pStore[2*i+1] )
            continue;
        uCanonPhase |= (1 << i);
        Temp = pStore[2*i+0];
        pStore[2*i+0] = pStore[2*i+1];
        pStore[2*i+1] = Temp;
        If_CluChangePhase( pIn, nVars, i );
    }
    // compute permutation
    for ( i = 0; i < nVars; i++ )
        pCanonPerm[i] = i;
    do {
        fChange = 0;
        for ( i = 0; i < nVars-1; i++ )
        {
            if ( pStore[2*i] <= pStore[2*(i+1)] )
                continue;
            Counter++;
            fChange = 1;

            Temp = pCanonPerm[i];
            pCanonPerm[i] = pCanonPerm[i+1];
            pCanonPerm[i+1] = Temp;

            Temp = pStore[2*i];
            pStore[2*i] = pStore[2*(i+1)];
            pStore[2*(i+1)] = Temp;

            Temp = pStore[2*i+1];
            pStore[2*i+1] = pStore[2*(i+1)+1];
            pStore[2*(i+1)+1] = Temp;

            If_CluSwapAdjacent( pOut, pIn, i, nVars );
            pTemp = pIn; pIn = pOut; pOut = pTemp;
        }
    } while ( fChange );
    // swap if it was moved an odd number of times
    if ( Counter & 1 )
        If_CluCopy( pOut, pIn, nVars );
    return uCanonPhase;
}
void If_CluSemiCanonicizeVerify( word * pTruth, word * pTruth0, int nVars, int * pCanonPerm, unsigned uCanonPhase )
{
    word pFunc[CLU_WRD_MAX], pGunc[CLU_WRD_MAX], * pIn = pTruth, * pOut = pFunc, * pTemp;
    int i, Temp, fChange, Counter = 0;
    If_CluCopy( pGunc, pTruth, nVars );
    // undo permutation
    do {
        fChange = 0;
        for ( i = 0; i < nVars-1; i++ )
        {
            if ( pCanonPerm[i] < pCanonPerm[i+1] )
                continue;

            Counter++;
            fChange = 1;

            Temp = pCanonPerm[i];
            pCanonPerm[i] = pCanonPerm[i+1];
            pCanonPerm[i+1] = Temp;

            If_CluSwapAdjacent( pOut, pIn, i, nVars );
            pTemp = pIn; pIn = pOut; pOut = pTemp;
        }
    } while ( fChange );
    if ( Counter & 1 )
        If_CluCopy( pOut, pIn, nVars );
    // undo phase
    for ( i = 0; i < nVars; i++ )
        if ( (uCanonPhase >> i) & 1 )
            If_CluChangePhase( pTruth, nVars, i );
    // compare
    if ( !If_CluEqual(pTruth0, pTruth, nVars) )
    {
        Kit_DsdPrintFromTruth( (unsigned*)pTruth0, nVars ); printf( "\n" );
        Kit_DsdPrintFromTruth( (unsigned*)pGunc, nVars ); printf( "\n" );
        Kit_DsdPrintFromTruth( (unsigned*)pTruth, nVars ); printf( "\n" );
        printf( "SemiCanonical verification FAILED!\n" );
    }
}


void If_CluPrintGroup( If_Grp_t * g )
{
    int i;
    printf( "Vars = %d   ", g->nVars );
    printf( "Myu = %d   {", g->nMyu );
    for ( i = 0; i < g->nVars; i++ )
        printf( " %c", 'a' + g->pVars[i] );
    printf( " }\n" );
}

void If_CluPrintConfig( int nVars, If_Grp_t * g, If_Grp_t * r, word BStruth, word * pFStruth )
{
    assert( r->nVars == nVars - g->nVars + 1 + (g->nMyu > 2) );
    If_CluPrintGroup( g );
    if ( g->nVars < 6 )
        BStruth = If_CluAdjust( BStruth, g->nVars );
    Kit_DsdPrintFromTruth( (unsigned *)&BStruth, g->nVars );
    printf( "\n" );
    If_CluPrintGroup( r );
    if ( r->nVars < 6 )
        pFStruth[0] = If_CluAdjust( pFStruth[0], r->nVars );
    Kit_DsdPrintFromTruth( (unsigned *)pFStruth, r->nVars );
    printf( "\n" );
}


void If_CluInitTruthTables()
{
    int i, k;
    assert( CLU_VAR_MAX <= 16 );
    for ( i = 0; i < 6; i++ )
        for ( k = 0; k < CLU_WRD_MAX; k++ )
            TruthAll[i][k] = Truth6[i];
    for ( i = 6; i < CLU_VAR_MAX; i++ )
        for ( k = 0; k < CLU_WRD_MAX; k++ )
            TruthAll[i][k] = ((k >> (i-6)) & 1) ? ~(word)0 : 0;

//    Extra_PrintHex( stdout, TruthAll[6], 8 ); printf( "\n" );
//    Extra_PrintHex( stdout, TruthAll[7], 8 ); printf( "\n" );
}


// verification
static void If_CluComposeLut( int nVars, If_Grp_t * g, word * t, word f[6][CLU_WRD_MAX], word * r )
{
    word c[CLU_WRD_MAX];
    int m, v;
    If_CluClear( r, nVars ); 
    for ( m = 0; m < (1<<g->nVars); m++ )
    {
        if ( !((t[m >> 6] >> (m & 63)) & 1) )
            continue;
        If_CluFill( c, nVars );
        for ( v = 0; v < g->nVars; v++ )
            if ( (m >> v) & 1 )
                If_CluAnd( c, c, f[v], nVars );
            else
                If_CluSharp( c, c, f[v], nVars );
        If_CluOr( r, r, c, nVars );
    }
}
void If_CluVerify( word * pF, int nVars, If_Grp_t * g, If_Grp_t * r, word BStruth, word * pFStruth )
{
    word pTTFans[6][CLU_WRD_MAX], pTTWire[CLU_WRD_MAX], pTTRes[CLU_WRD_MAX];
    int i;
    assert( g->nVars <= 6 && r->nVars <= 6 );

    if ( TruthAll[0][0] == 0 )
        If_CluInitTruthTables();

    for ( i = 0; i < g->nVars; i++ )
        If_CluCopy( pTTFans[i], TruthAll[(int)g->pVars[i]], nVars );
    If_CluComposeLut( nVars, g, &BStruth, pTTFans, pTTWire );

    for ( i = 0; i < r->nVars; i++ )
        if ( r->pVars[i] == nVars )
            If_CluCopy( pTTFans[i], pTTWire, nVars );
        else
            If_CluCopy( pTTFans[i], TruthAll[(int)r->pVars[i]], nVars );
    If_CluComposeLut( nVars, r, pFStruth, pTTFans, pTTRes );

    if ( !If_CluEqual(pTTRes, pF, nVars) )
    {
        printf( "\n" );
        If_CluPrintConfig( nVars, g, r, BStruth, pFStruth );
        Kit_DsdPrintFromTruth( (unsigned*)pTTRes, nVars ); printf( "\n" );
        Kit_DsdPrintFromTruth( (unsigned*)pF, nVars ); printf( "\n" );
//        Extra_PrintHex( stdout, (unsigned *)pF, nVars ); printf( "\n" );
        printf( "Verification FAILED!\n" );
    }
    else if ( 0 )
    {
        printf( "\n" );
        If_CluPrintConfig( nVars, g, r, BStruth, pFStruth );
        //Kit_DsdPrintFromTruth( (unsigned*)pTTRes, nVars ); printf( "\n" );
        //Kit_DsdPrintFromTruth( (unsigned*)pF, nVars ); printf( "\n" );        
        printf( "Verification succeed!\n" );
    }
}
void If_CluVerify3( word * pF, int nVars, If_Grp_t * g, If_Grp_t * g2, If_Grp_t * r, word BStruth, word BStruth2, word FStruth )
{
    word pTTFans[6][CLU_WRD_MAX], pTTWire[CLU_WRD_MAX], pTTWire2[CLU_WRD_MAX], pTTRes[CLU_WRD_MAX];
    int i;
    assert( g->nVars >= 2 && g2->nVars >= 2 && r->nVars >= 2 );
    assert( g->nVars <= 6 && g2->nVars <= 6 && r->nVars <= 6 );

    if ( TruthAll[0][0] == 0 )
        If_CluInitTruthTables();

    for ( i = 0; i < g->nVars; i++ )
        If_CluCopy( pTTFans[i], TruthAll[(int)g->pVars[i]], nVars );
    If_CluComposeLut( nVars, g, &BStruth, pTTFans, pTTWire );

    for ( i = 0; i < g2->nVars; i++ )
        If_CluCopy( pTTFans[i], TruthAll[(int)g2->pVars[i]], nVars );
    If_CluComposeLut( nVars, g2, &BStruth2, pTTFans, pTTWire2 );

    for ( i = 0; i < r->nVars; i++ )
        if ( r->pVars[i] == nVars )
            If_CluCopy( pTTFans[i], pTTWire, nVars );
        else if ( r->pVars[i] == nVars + 1 )
            If_CluCopy( pTTFans[i], pTTWire2, nVars );
        else
            If_CluCopy( pTTFans[i], TruthAll[(int)r->pVars[i]], nVars );
    If_CluComposeLut( nVars, r, &FStruth, pTTFans, pTTRes );

    if ( !If_CluEqual(pTTRes, pF, nVars) )
    {
        printf( "%d\n", nVars );
//        If_CluPrintConfig( nVars, g, r, BStruth, pFStruth );
//        Extra_PrintHex( stdout, (unsigned *)pF, nVars ); printf( "\n" );

        Kit_DsdPrintFromTruth( (unsigned*)&BStruth, g->nVars );   printf( "    " ); If_CluPrintGroup(g);  // printf( "\n" );
        Kit_DsdPrintFromTruth( (unsigned*)&BStruth2, g2->nVars ); printf( "    " ); If_CluPrintGroup(g2); // printf( "\n" );
        Kit_DsdPrintFromTruth( (unsigned*)&FStruth, r->nVars );   printf( "    " ); If_CluPrintGroup(r);  // printf( "\n" );

        Kit_DsdPrintFromTruth( (unsigned*)pTTWire, nVars ); printf( "\n" );
        Kit_DsdPrintFromTruth( (unsigned*)pTTWire2, nVars ); printf( "\n" );
        Kit_DsdPrintFromTruth( (unsigned*)pTTRes, nVars ); printf( "\n" );
        Kit_DsdPrintFromTruth( (unsigned*)pF, nVars ); printf( "\n" );
//        Extra_PrintHex( stdout, (unsigned *)pF, nVars ); printf( "\n" );
        printf( "Verification FAILED!\n" );
    }
//    else
//        printf( "Verification succeed!\n" );
}



void If_CluSwapVars( word * pTruth, int nVars, int * V2P, int * P2V, int iVar, int jVar )
{
    word low2High, high2Low, temp;
    int nWords = If_CluWordNum(nVars);
    int shift, step, iStep, jStep;
    int w = 0, i = 0, j = 0;
    static word PPMasks[6][6] = {
        { ABC_CONST(0x2222222222222222), ABC_CONST(0x0A0A0A0A0A0A0A0A), ABC_CONST(0x00AA00AA00AA00AA), ABC_CONST(0x0000AAAA0000AAAA), ABC_CONST(0x00000000AAAAAAAA), ABC_CONST(0xAAAAAAAAAAAAAAAA) },
        { ABC_CONST(0x0000000000000000), ABC_CONST(0x0C0C0C0C0C0C0C0C), ABC_CONST(0x00CC00CC00CC00CC), ABC_CONST(0x0000CCCC0000CCCC), ABC_CONST(0x00000000CCCCCCCC), ABC_CONST(0xCCCCCCCCCCCCCCCC) },
        { ABC_CONST(0x0000000000000000), ABC_CONST(0x0000000000000000), ABC_CONST(0x00F000F000F000F0), ABC_CONST(0x0000F0F00000F0F0), ABC_CONST(0x00000000F0F0F0F0), ABC_CONST(0xF0F0F0F0F0F0F0F0) },
        { ABC_CONST(0x0000000000000000), ABC_CONST(0x0000000000000000), ABC_CONST(0x0000000000000000), ABC_CONST(0x0000FF000000FF00), ABC_CONST(0x00000000FF00FF00), ABC_CONST(0xFF00FF00FF00FF00) },
        { ABC_CONST(0x0000000000000000), ABC_CONST(0x0000000000000000), ABC_CONST(0x0000000000000000), ABC_CONST(0x0000000000000000), ABC_CONST(0x00000000FFFF0000), ABC_CONST(0xFFFF0000FFFF0000) },
        { ABC_CONST(0x0000000000000000), ABC_CONST(0x0000000000000000), ABC_CONST(0x0000000000000000), ABC_CONST(0x0000000000000000), ABC_CONST(0x0000000000000000), ABC_CONST(0xFFFFFFFF00000000) }
    };
    if( iVar == jVar )
        return;
    if( jVar < iVar )
    {
        int varTemp = jVar;
        jVar = iVar;
        iVar = varTemp;
    }
    if ( iVar <= 5 && jVar <= 5 )
    {
        shift = (1 <<  jVar) - (1 << iVar);
        for ( w = 0; w < nWords; w++ )
        {
            low2High = (pTruth[w] & PPMasks[iVar][jVar - 1] ) << shift;
            pTruth[w] &= ~PPMasks[iVar][jVar - 1];
            high2Low = (pTruth[w] & (PPMasks[iVar][jVar - 1] << shift )) >> shift;
            pTruth[w] &= ~ (PPMasks[iVar][jVar - 1] << shift);
            pTruth[w] = pTruth[w] | low2High | high2Low;
        }
    }
    else if( iVar <= 5 && jVar > 5 )
    {
        step = If_CluWordNum(jVar + 1)/2;
        shift = 1 << iVar;
        for ( w = 0; w < nWords; w += 2*step )
        {
            for (j = 0; j < step; j++)
            {
                low2High = (pTruth[w + j] & PPMasks[iVar][5]) >> shift;
                pTruth[w + j] &= ~PPMasks[iVar][5];
                high2Low = (pTruth[w + step + j] & (PPMasks[iVar][5] >> shift)) << shift;
                pTruth[w + step + j] &= ~(PPMasks[iVar][5] >> shift);
                pTruth[w + j] |= high2Low;
                pTruth[w + step + j] |= low2High;            
            }
        }
    }
    else
    {
        iStep = If_CluWordNum(iVar + 1)/2;
        jStep = If_CluWordNum(jVar + 1)/2;
        for (w = 0; w < nWords; w += 2*jStep)
        {
            for (i = 0; i < jStep; i += 2*iStep)
            {
                for (j = 0; j < iStep; j++)
                {
                    temp = pTruth[w + iStep + i + j];
                    pTruth[w + iStep + i + j] = pTruth[w + jStep + i + j];
                    pTruth[w + jStep + i + j] = temp;
                }
            }
        }
    }    
    if ( V2P && P2V )
    {
        V2P[P2V[iVar]] = jVar;
        V2P[P2V[jVar]] = iVar;
        P2V[iVar] ^= P2V[jVar];
        P2V[jVar] ^= P2V[iVar];
        P2V[iVar] ^= P2V[jVar];
    }
}
void If_CluReverseOrder( word * pTruth, int nVars, int * V2P, int * P2V, int iVarStart )
{
    int i, j, k;
    for ( k = 0; k < (nVars-iVarStart)/2 ; k++ )
    {
        i = iVarStart + k;
        j = nVars - 1 - k;
        If_CluSwapVars( pTruth, nVars, V2P, P2V, i, j );
    }
}

// moves one var (v) to the given position (p)
void If_CluMoveVar2( word * pF, int nVars, int * Var2Pla, int * Pla2Var, int v, int p )
{
    If_CluSwapVars( pF, nVars, Var2Pla, Pla2Var, Var2Pla[v], p );
}

// moves one var (v) to the given position (p)
void If_CluMoveVar( word * pF, int nVars, int * Var2Pla, int * Pla2Var, int v, int p )
{
    word pG[CLU_WRD_MAX], * pIn = pF, * pOut = pG, * pTemp;
    int iPlace0, iPlace1, Count = 0;
    assert( v >= 0 && v < nVars );
    while ( Var2Pla[v] < p )
    {
        iPlace0 = Var2Pla[v];
        iPlace1 = Var2Pla[v]+1;
        If_CluSwapAdjacent( pOut, pIn, iPlace0, nVars );
        pTemp = pIn; pIn = pOut, pOut = pTemp;
        Var2Pla[Pla2Var[iPlace0]]++;
        Var2Pla[Pla2Var[iPlace1]]--;
        Pla2Var[iPlace0] ^= Pla2Var[iPlace1];
        Pla2Var[iPlace1] ^= Pla2Var[iPlace0];
        Pla2Var[iPlace0] ^= Pla2Var[iPlace1];
        Count++;
    }
    while ( Var2Pla[v] > p )
    {
        iPlace0 = Var2Pla[v]-1;
        iPlace1 = Var2Pla[v];
        If_CluSwapAdjacent( pOut, pIn, iPlace0, nVars );
        pTemp = pIn; pIn = pOut, pOut = pTemp;
        Var2Pla[Pla2Var[iPlace0]]++;
        Var2Pla[Pla2Var[iPlace1]]--;
        Pla2Var[iPlace0] ^= Pla2Var[iPlace1];
        Pla2Var[iPlace1] ^= Pla2Var[iPlace0];
        Pla2Var[iPlace0] ^= Pla2Var[iPlace1];
        Count++;
    }
    if ( Count & 1 )
        If_CluCopy( pF, pIn, nVars );
    assert( Pla2Var[p] == v );
}

// moves vars to be the most signiticant ones (Group[0] is MSB)
void If_CluMoveGroupToMsb( word * pF, int nVars, int * V2P, int * P2V, If_Grp_t * g )
{
    int v;
    for ( v = 0; v < g->nVars; v++ )
        If_CluMoveVar( pF, nVars, V2P, P2V, g->pVars[g->nVars - 1 - v], nVars - 1 - v );
}


// reverses the variable order
void If_CluReverseOrder_old( word * pF, int nVars, int * V2P, int * P2V, int iVarStart )
{
    word pG[CLU_WRD_MAX];
    int v;

    If_CluCopy( pG, pF, nVars );

//    for ( v = 0; v < nVars; v++ )
//        printf( "%c ", 'a' + P2V[v] );
//    printf( "  ---  " );

    for ( v = iVarStart; v < nVars; v++ )
        If_CluMoveVar( pF, nVars, V2P, P2V, P2V[iVarStart], nVars - 1 - (v - iVarStart) );

//    for ( v = 0; v < nVars; v++ )
//        printf( "%c ", 'a' + P2V[v] );
//    printf( "\n" );

//    if ( iVarStart > 0 )
//        return;

    If_CluReverseOrder( pG, nVars, NULL, NULL, iVarStart );
    if ( If_CluEqual( pG, pF, nVars ) )
    {
//        printf( "+" );
    }
    else
    {
/*
        printf( "\n" );
        Kit_DsdPrintFromTruth( (unsigned*)pF, nVars ); printf( "\n" );
        Kit_DsdPrintFromTruth( (unsigned*)pG, nVars ); 
        printf( "\n" );
*/
        printf( "%d ", nVars );
    }
}

// return the number of cofactors w.r.t. the topmost vars (nBSsize)
int If_CluCountCofs( word * pF, int nVars, int nBSsize, int iShift, word pCofs[3][CLU_WRD_MAX/4] )
{
    word iCofs[128] = {0}, iCof, Result = 0;
    word * pCofA, * pCofB;
    int nMints = (1 << nBSsize);
    int i, c, w, nCofs;
    assert( nBSsize >= 2 && nBSsize <= 7 && nBSsize < nVars );
    if ( nVars - nBSsize < 6 )
    {
        int nShift = (1 << (nVars - nBSsize));
        word Mask  = ((((word)1) << nShift) - 1);
        for ( nCofs = i = 0; i < nMints; i++ )
        {
            iCof = (pF[(iShift + i * nShift) / 64] >> ((iShift + i * nShift) & 63)) & Mask;
            for ( c = 0; c < nCofs; c++ )
                if ( iCof == iCofs[c] )
                    break;
            if ( c == nCofs )
                iCofs[nCofs++] = iCof;
            if ( pCofs && iCof != iCofs[0] )
                Result |= (((word)1) << i);
            if ( nCofs == 5 )
                break;
        }
        if ( nCofs <= 2 && pCofs )
        {
            assert( nBSsize <= 6 );
            pCofs[0][0] = iCofs[0];
            pCofs[1][0] = (nCofs == 2) ? iCofs[1] : iCofs[0];
            pCofs[2][0] = Result;
        }
    }
    else
    {
        int nWords = If_CluWordNum( nVars - nBSsize );
        assert( nWords * nMints == If_CluWordNum(nVars) );
        for ( nCofs = i = 0; i < nMints; i++ )
        {
            pCofA = pF + i * nWords;
            for ( c = 0; c < nCofs; c++ )
            {
                pCofB = pF + iCofs[c] * nWords;
                for ( w = 0; w < nWords; w++ )
                    if ( pCofA[w] != pCofB[w] )
                        break;
                if ( w == nWords )
                    break;
            }
            if ( c == nCofs )
                iCofs[nCofs++] = i;
            if ( pCofs )
            {
                assert( nBSsize <= 6 );
                pCofB = pF + iCofs[0] * nWords;
                for ( w = 0; w < nWords; w++ )
                    if ( pCofA[w] != pCofB[w] )
                        break;
                if ( w != nWords )
                    Result |= (((word)1) << i);
            }
            if ( nCofs == 5 )
                break;
        }
        if ( nCofs <= 2 && pCofs )
        {
            If_CluCopy( pCofs[0], pF + iCofs[0] * nWords, nVars - nBSsize );
            If_CluCopy( pCofs[1], pF + ((nCofs == 2) ? iCofs[1] : iCofs[0]) * nWords, nVars - nBSsize );
            pCofs[2][0] = Result;
        }
    }
    assert( nCofs >= 1 && nCofs <= 5 );
    return nCofs;
}

// return the number of cofactors w.r.t. the topmost vars (nBSsize)
int If_CluCountCofs4( word * pF, int nVars, int nBSsize, word pCofs[6][CLU_WRD_MAX/4] )
{
    word iCofs[128] = {0}, iCof, Result0 = 0, Result1 = 0;
    int nMints = (1 << nBSsize);
    int i, c, nCofs = 0;
    assert( pCofs );
    assert( nBSsize >= 2 && nBSsize <= 6 && nBSsize < nVars );
    if ( nVars - nBSsize < 6 )
    {
        int nShift = (1 << (nVars - nBSsize));
        word Mask  = ((((word)1) << nShift) - 1);
        for ( nCofs = i = 0; i < nMints; i++ )
        {
            iCof = (pF[(i * nShift) / 64] >> ((i * nShift) & 63)) & Mask;
            for ( c = 0; c < nCofs; c++ )
                if ( iCof == iCofs[c] )
                    break;
            if ( c == nCofs )
                iCofs[nCofs++] = iCof;
            if ( c == 1 || c == 3 )
                Result0 |= (((word)1) << i);
            if ( c == 2 || c == 3 )
                Result1 |= (((word)1) << i);
        }
        assert( nCofs >= 3 && nCofs <= 4 );
        pCofs[0][0] = iCofs[0];
        pCofs[1][0] = iCofs[1];
        pCofs[2][0] = iCofs[2];
        pCofs[3][0] = (nCofs == 4) ? iCofs[3] : iCofs[2];
        pCofs[4][0] = Result0;
        pCofs[5][0] = Result1;
    }
    else
        assert( 0 );
    return nCofs;
}

void If_CluCofactors( word * pF, int nVars, int iVar, word * pCof0, word * pCof1 )
{
    int nWords = If_CluWordNum( nVars );
    assert( iVar < nVars );
    if ( iVar < 6 )
    {
        int i, Shift = (1 << iVar);
        for ( i = 0; i < nWords; i++ )
        {
            pCof0[i] = (pF[i] & ~Truth6[iVar]) | ((pF[i] & ~Truth6[iVar]) << Shift);
            pCof1[i] = (pF[i] &  Truth6[iVar]) | ((pF[i] &  Truth6[iVar]) >> Shift);
        }
    }
    else
    {
        int i, k, Step = (1 << (iVar - 6));
        for ( k = 0; k < nWords; k += 2*Step )
        {
            for ( i = 0; i < Step; i++ )
            {
                pCof0[i] = pCof0[Step+i] = pF[i];
                pCof1[i] = pCof1[Step+i] = pF[Step+i];
            }
            pF    += 2*Step;
            pCof0 += 2*Step;
            pCof1 += 2*Step;
        }
    }
}

// returns 1 if we have special case of cofactors; otherwise, returns 0
int If_CluDetectSpecialCaseCofs( word * pF, int nVars, int iVar )
{
    word Cof0, Cof1;
    int State[6] = {0};
    int i, nWords = If_CluWordNum( nVars );
    assert( iVar < nVars );
    if ( iVar < 6 )
    {
        int Shift = (1 << iVar);
        for ( i = 0; i < nWords; i++ )
        {
            Cof0 =  (pF[i] & ~Truth6[iVar]);
            Cof1 = ((pF[i] &  Truth6[iVar]) >> Shift);

            if ( Cof0 == 0 )
                State[0]++;
            else if ( Cof0 == ~Truth6[iVar] )
                State[1]++;
            else if ( Cof1 == 0 )
                State[2]++;
            else if ( Cof1 == ~Truth6[iVar] )
                State[3]++;
            else if ( Cof0 == ~Cof1 )
                State[4]++;
            else if ( Cof0 == Cof1 )
                State[5]++;
        }
    }
    else
    {
        int k, Step = (1 << (iVar - 6));
        for ( k = 0; k < nWords; k += 2*Step )
        {
            for ( i = 0; i < Step; i++ )
            {
                Cof0 = pF[i];
                Cof1 = pF[Step+i];

                if ( Cof0 == 0 )
                    State[0]++;
                else if ( Cof0 == ~(word)0 )
                    State[1]++;
                else if ( Cof1 == 0 )
                    State[2]++;
                else if ( Cof1 == ~(word)0 )
                    State[3]++;
                else if ( Cof0 == ~Cof1 )
                    State[4]++;
                else if ( Cof0 == Cof1 )
                    State[5]++;
            }
            pF    += 2*Step;
        }
        nWords /= 2;
    }
    assert( State[5] != nWords );
    for ( i = 0; i < 5; i++ )
    {
        assert( State[i] <= nWords );
        if ( State[i] == nWords )
            return i;
    }
    return -1;
}

// returns 1 if we have special case of cofactors; otherwise, returns 0
If_Grp_t If_CluDecUsingCofs( word * pTruth, int nVars, int nLutLeaf )
{
    If_Grp_t G = {0};
    word pF2[CLU_WRD_MAX], * pF = pF2;
    int Var2Pla[CLU_VAR_MAX+2], Pla2Var[CLU_VAR_MAX+2];
    int V2P[CLU_VAR_MAX+2], P2V[CLU_VAR_MAX+2];
    int nVarsNeeded = nVars - nLutLeaf;
    int v, i, k, iVar, State;
//Kit_DsdPrintFromTruth( (unsigned*)pTruth, nVars ); printf( "\n" );
    // create local copy
    If_CluCopy( pF, pTruth, nVars );
    for ( k = 0; k < nVars; k++ )
        Var2Pla[k] = Pla2Var[k] = k;
    // find decomposable vars 
    for ( i = 0; i < nVarsNeeded; i++ )
    {
        for ( v = nVars - 1; v >= 0; v-- )
        {
            State = If_CluDetectSpecialCaseCofs( pF, nVars, v );
            if ( State == -1 )
                continue;
            // update the variable place
            iVar = Pla2Var[v];
            while ( Var2Pla[iVar] < nVars - 1 )
            {
                int iPlace0 = Var2Pla[iVar];
                int iPlace1 = Var2Pla[iVar]+1;
                Var2Pla[Pla2Var[iPlace0]]++;
                Var2Pla[Pla2Var[iPlace1]]--;
                Pla2Var[iPlace0] ^= Pla2Var[iPlace1];
                Pla2Var[iPlace1] ^= Pla2Var[iPlace0];
                Pla2Var[iPlace0] ^= Pla2Var[iPlace1];
            }
            // move this variable to the top
            for ( k = 0; k < nVars; k++ )
                V2P[k] = P2V[k] = k;
//Kit_DsdPrintFromTruth( (unsigned*)pF, nVars ); printf( "\n" );
            If_CluMoveVar( pF, nVars, V2P, P2V, v, nVars - 1 );
//Kit_DsdPrintFromTruth( (unsigned*)pF, nVars ); printf( "\n" );
            // choose cofactor to follow
            iVar = nVars - 1;
            if ( State == 0 || State == 1 ) // need cof1
            {
                if ( iVar < 6 )
                    pF[0] = (pF[0] &  Truth6[iVar]) | ((pF[0] &  Truth6[iVar]) >> (1 << iVar));
                else
                    pF += If_CluWordNum( nVars ) / 2;
            }
            else // need cof0
            {
                if ( iVar < 6 )
                    pF[0] = (pF[0] & ~Truth6[iVar]) | ((pF[0] & ~Truth6[iVar]) << (1 << iVar));
            }
            // update the variable count
            nVars--;
            break;
        }
        if ( v == -1 )
            return G;
    }
    // create the resulting group
    G.nVars = nLutLeaf;
    G.nMyu = 2;
    for ( v = 0; v < G.nVars; v++ )
        G.pVars[v] = Pla2Var[v];
    return G;
}



// deriving decomposition
word If_CluDeriveDisjoint( word * pF, int nVars, int * V2P, int * P2V, If_Grp_t * g, If_Grp_t * r )
{
    word pCofs[3][CLU_WRD_MAX/4];
    int i, RetValue, nFSset = nVars - g->nVars;
    RetValue = If_CluCountCofs( pF, nVars, g->nVars, 0, pCofs );
//    assert( RetValue == 2 );

    if ( nFSset < 6 )
        pF[0] = (pCofs[1][0] << (1 << nFSset)) | pCofs[0][0];
    else
    {
        If_CluCopy( pF, pCofs[0], nFSset );
        If_CluCopy( pF + If_CluWordNum(nFSset), pCofs[1], nFSset );
    }
    // create the resulting group
    if ( r )
    {
        r->nVars = nFSset + 1;
        r->nMyu = 0;
        for ( i = 0; i < nFSset; i++ )
            r->pVars[i] = P2V[i];
        r->pVars[nFSset] = nVars;
    }
    return pCofs[2][0];
}
void If_CluDeriveDisjoint4( word * pF, int nVars, int * V2P, int * P2V, If_Grp_t * g, If_Grp_t * r, word * pTruth0, word * pTruth1 )
{
    word pCofs[6][CLU_WRD_MAX/4];
    word Cof0, Cof1;
    int i, RetValue, nFSset = nVars - g->nVars;

    assert( g->nVars <= 6 && nFSset <= 4 );

    RetValue = If_CluCountCofs4( pF, nVars, g->nVars, pCofs );
    if ( RetValue != 3 && RetValue != 4 )
        printf( "If_CluDeriveDisjoint4(): Error!!!\n" );

    Cof0  = (pCofs[1][0] << (1 << nFSset)) | pCofs[0][0];
    Cof1  = (pCofs[3][0] << (1 << nFSset)) | pCofs[2][0];
    pF[0] = (Cof1 << (1 << (nFSset+1))) | Cof0;
    pF[0] = If_CluAdjust( pF[0], nFSset + 2 );

    // create the resulting group
    r->nVars = nFSset + 2;
    r->nMyu = 0;
    for ( i = 0; i < nFSset; i++ )
        r->pVars[i] = P2V[i];
    r->pVars[nFSset] = nVars;
    r->pVars[nFSset+1] = nVars+1;

    *pTruth0 = If_CluAdjust( pCofs[4][0], g->nVars );
    *pTruth1 = If_CluAdjust( pCofs[5][0], g->nVars );
}

word If_CluDeriveNonDisjoint( word * pF, int nVars, int * V2P, int * P2V, If_Grp_t * g, If_Grp_t * r )
{
    word pCofs[2][CLU_WRD_MAX];
    word Truth0, Truth1, Truth;
    int i, nFSset = nVars - g->nVars, nFSset1 = nFSset + 1;
    If_CluCofactors( pF, nVars, nVars - 1, pCofs[0], pCofs[1] );

//    Extra_PrintHex( stdout, (unsigned *)pCofs[0], nVars ); printf( "\n" );
//    Extra_PrintHex( stdout, (unsigned *)pCofs[1], nVars ); printf( "\n" );

    g->nVars--;
    Truth0 = If_CluDeriveDisjoint( pCofs[0], nVars - 1, V2P, P2V, g, NULL );
    Truth1 = If_CluDeriveDisjoint( pCofs[1], nVars - 1, V2P, P2V, g, NULL );
    Truth  = (Truth1 << (1 << g->nVars)) | Truth0;
    g->nVars++;
    if ( nFSset1 < 6 )
        pF[0] = (pCofs[1][0] << (1 << nFSset1)) | pCofs[0][0];
    else
    {
        If_CluCopy( pF, pCofs[0], nFSset1 );
        If_CluCopy( pF + If_CluWordNum(nFSset1), pCofs[1], nFSset1 );
    }

//    Extra_PrintHex( stdout, (unsigned *)&Truth0, 6 ); printf( "\n" );
//    Extra_PrintHex( stdout, (unsigned *)&Truth1, 6 ); printf( "\n" );
//    Extra_PrintHex( stdout, (unsigned *)&pCofs[0][0], 6 ); printf( "\n" );
//    Extra_PrintHex( stdout, (unsigned *)&pCofs[1][0], 6 ); printf( "\n" );
//    Extra_PrintHex( stdout, (unsigned *)&Truth, 6 ); printf( "\n" );
//    Extra_PrintHex( stdout, (unsigned *)&pF[0], 6 ); printf( "\n" );

    // create the resulting group
    r->nVars = nFSset + 2;
    r->nMyu = 0;
    for ( i = 0; i < nFSset; i++ )
        r->pVars[i] = P2V[i];
    r->pVars[nFSset] = nVars;
    r->pVars[nFSset+1] = g->pVars[g->nVars - 1];
    return Truth;
}

// check non-disjoint decomposition
int If_CluCheckNonDisjointGroup( word * pF, int nVars, int * V2P, int * P2V, If_Grp_t * g )
{
    int v, i, nCofsBest2;
    if ( (g->nMyu == 3 || g->nMyu == 4) )
    {
        word pCofs[2][CLU_WRD_MAX];
        // try cofactoring w.r.t. each variable
        for ( v = 0; v < g->nVars; v++ )
        {
            If_CluCofactors( pF, nVars, V2P[(int)g->pVars[v]], pCofs[0], pCofs[1] );
            nCofsBest2 = If_CluCountCofs( pCofs[0], nVars, g->nVars, 0, NULL );
            if ( nCofsBest2 > 2 )
                continue;
            nCofsBest2 = If_CluCountCofs( pCofs[1], nVars, g->nVars, 0, NULL );
            if ( nCofsBest2 > 2 )
                continue;
            // found good shared variable - move to the end
            If_CluMoveVar( pF, nVars, V2P, P2V, g->pVars[v], nVars-1 );
            for ( i = 0; i < g->nVars; i++ )
                g->pVars[i] = P2V[nVars-g->nVars+i];
            return 1;
        }
    }
    return 0;
}


// finds a good var group (cof count < 6; vars are MSBs)
If_Grp_t If_CluFindGroup( word * pF, int nVars, int iVarStart, int iVarStop, int * V2P, int * P2V, int nBSsize, int fDisjoint )
{
    int fVerbose = 0;
    int nRounds = 2;//nBSsize;
    If_Grp_t G = {0}, * g = &G;//, BestG = {0};
    int i, r, v, nCofs, VarBest, nCofsBest2;
    assert( nVars > nBSsize && nVars >= nBSsize + iVarStart && nVars <= CLU_VAR_MAX );
    assert( nBSsize >= 2 && nBSsize <= 6 );
    assert( !iVarStart || !iVarStop );
    // start with the default group
    g->nVars = nBSsize;
    g->nMyu = If_CluCountCofs( pF, nVars, nBSsize, 0, NULL );
    for ( i = 0; i < nBSsize; i++ )
        g->pVars[i] = P2V[nVars-nBSsize+i];
    // check if good enough
    if ( g->nMyu == 2 )
        return G;
    if ( !fDisjoint && If_CluCheckNonDisjointGroup( pF, nVars, V2P, P2V, g ) )
    {
//        BestG = G;
        return G;
    }
    if ( nVars == nBSsize + iVarStart )
    {
        g->nVars = 0;
        return G;
    }

    if ( fVerbose )
    {
        printf( "Iter %2d  ", -1 );
        If_CluPrintGroup( g );
    }

    // try to find better group
    for ( r = 0; r < nRounds; r++ )
    {
        if ( nBSsize < nVars-1 )
        {
            // find the best var to add
            VarBest = P2V[nVars-1-nBSsize];
            nCofsBest2 = If_CluCountCofs( pF, nVars, nBSsize+1, 0, NULL );
            for ( v = nVars-2-nBSsize; v >= iVarStart; v-- )
            {
//                If_CluMoveVar( pF, nVars, V2P, P2V, P2V[v], nVars-1-nBSsize );
                If_CluMoveVar2( pF, nVars, V2P, P2V, P2V[v], nVars-1-nBSsize );
                nCofs = If_CluCountCofs( pF, nVars, nBSsize+1, 0, NULL );
                if ( nCofsBest2 >= nCofs )
                {
                    nCofsBest2 = nCofs;
                    VarBest = P2V[nVars-1-nBSsize];
                }
            }
            // go back
//            If_CluMoveVar( pF, nVars, V2P, P2V, VarBest, nVars-1-nBSsize );
            If_CluMoveVar2( pF, nVars, V2P, P2V, VarBest, nVars-1-nBSsize );
            // update best bound set
            nCofs = If_CluCountCofs( pF, nVars, nBSsize+1, 0, NULL );
            assert( nCofs == nCofsBest2 );
        }

        // find the best var to remove
        VarBest = P2V[nVars-1-nBSsize];
        nCofsBest2 = If_CluCountCofs( pF, nVars, nBSsize, 0, NULL );
        for ( v = nVars-nBSsize; v < nVars-iVarStop; v++ )
        {
//            If_CluMoveVar( pF, nVars, V2P, P2V, P2V[v], nVars-1-nBSsize );
            If_CluMoveVar2( pF, nVars, V2P, P2V, P2V[v], nVars-1-nBSsize );
            nCofs = If_CluCountCofs( pF, nVars, nBSsize, 0, NULL );
            if ( nCofsBest2 >= nCofs )
            {
                nCofsBest2 = nCofs;
                VarBest = P2V[nVars-1-nBSsize];
            }
        }

        // go back
//        If_CluMoveVar( pF, nVars, V2P, P2V, VarBest, nVars-1-nBSsize );
        If_CluMoveVar2( pF, nVars, V2P, P2V, VarBest, nVars-1-nBSsize );
        // update best bound set
        nCofs = If_CluCountCofs( pF, nVars, nBSsize, 0, NULL );
        assert( nCofs == nCofsBest2 );
        if ( g->nMyu >= nCofs )
        {
            g->nVars = nBSsize;
            g->nMyu = nCofs;
            for ( i = 0; i < nBSsize; i++ )
                g->pVars[i] = P2V[nVars-nBSsize+i];
        }

        if ( fVerbose )
        {
            printf( "Iter %2d  ", r );
            If_CluPrintGroup( g );
        }

        // check if good enough
        if ( g->nMyu == 2 )
            return G;
        if ( !fDisjoint && If_CluCheckNonDisjointGroup( pF, nVars, V2P, P2V, g ) )
        {
//            BestG = G;
            return G;
        }
    }

    assert( r == nRounds );
    g->nVars = 0;
    return G;
//    return BestG;
}


// double check that the given group has a decomposition
void If_CluCheckGroup( word * pTruth, int nVars, If_Grp_t * g )
{
    word pF[CLU_WRD_MAX];
    int v, nCofs, V2P[CLU_VAR_MAX], P2V[CLU_VAR_MAX];
    assert( g->nVars >= 2 && g->nVars <= 6 ); // vars
    assert( g->nMyu >= 2 && g->nMyu <= 4 ); // cofs
    // create permutation
    for ( v = 0; v < nVars; v++ )
        V2P[v] = P2V[v] = v;
    // create truth table
    If_CluCopy( pF, pTruth, nVars );
    // move group up
    If_CluMoveGroupToMsb( pF, nVars, V2P, P2V, g );
    // check the number of cofactors
    nCofs = If_CluCountCofs( pF, nVars, g->nVars, 0, NULL );
    if ( nCofs != g->nMyu )
        printf( "Group check 0 has failed.\n" );
    // check compatible
    if ( nCofs > 2 )
    {
        nCofs = If_CluCountCofs( pF, nVars-1, g->nVars-1, 0, NULL );
        if ( nCofs > 2 )
            printf( "Group check 1 has failed.\n" );
        nCofs = If_CluCountCofs( pF, nVars-1, g->nVars-1, (1 << (nVars-1)), NULL );
        if ( nCofs > 2 )
            printf( "Group check 2 has failed.\n" );
    }
}


// double check that the permutation derived is correct
void If_CluCheckPerm( word * pTruth, word * pF, int nVars, int * V2P, int * P2V )
{
    int i;
    for ( i = 0; i < nVars; i++ )
        If_CluMoveVar( pF, nVars, V2P, P2V, i, i );

    if ( !If_CluEqual( pTruth, pF, nVars ) )
        printf( "Permutation FAILED.\n" );
//    else
//        printf( "Permutation successful\n" );
}




static inline int If_CluSuppIsMinBase( int Supp )
{
    return (Supp & (Supp+1)) == 0;
}
static inline int If_CluHasVar( word * t, int nVars, int iVar )
{
    int nWords = If_CluWordNum( nVars );
    assert( iVar < nVars );
    if ( iVar < 6 )
    {
        int i, Shift = (1 << iVar);
        for ( i = 0; i < nWords; i++ )
            if ( (t[i] & ~Truth6[iVar]) != ((t[i] & Truth6[iVar]) >> Shift) )
                return 1;
        return 0;
    }
    else
    {
        int i, k, Step = (1 << (iVar - 6));
        for ( k = 0; k < nWords; k += 2*Step )
        {
            for ( i = 0; i < Step; i++ )
                if ( t[i] != t[Step+i] )
                    return 1;
            t += 2*Step;
        }
        return 0;
    }
}
static inline int If_CluSupport( word * t, int nVars )
{
    int v, Supp = 0;
    for ( v = 0; v < nVars; v++ )
        if ( If_CluHasVar( t, nVars, v ) )
            Supp |= (1 << v);
    return Supp;
}
int If_CluSupportSize( word * t, int nVars )
{
    int v, SuppSize = 0;
    for ( v = 0; v < nVars; v++ )
        if ( If_CluHasVar( t, nVars, v ) )
            SuppSize++;
    return SuppSize;
}
static inline void If_CluTruthShrink( word * pF, int nVars, int nVarsAll, unsigned Phase )
{
    word pG[CLU_WRD_MAX], * pIn = pF, * pOut = pG, * pTemp;
    int i, k, Var = 0, Counter = 0;
    assert( nVarsAll <= 16 );
    for ( i = 0; i < nVarsAll; i++ )
        if ( Phase & (1 << i) )
        {
            for ( k = i-1; k >= Var; k-- )
            {
                If_CluSwapAdjacent( pOut, pIn, k, nVarsAll );
                pTemp = pIn; pIn = pOut, pOut = pTemp;
                Counter++;
            }
            Var++;
        }
    assert( Var == nVars );
    // swap if it was moved an odd number of times
    if ( Counter & 1 )
        If_CluCopy( pOut, pIn, nVarsAll );
}
int If_CluMinimumBase( word * t, int * pSupp, int nVarsAll, int * pnVars )
{
    int v, iVar = 0, uSupp = 0;
    assert( nVarsAll <= 16 );
    for ( v = 0; v < nVarsAll; v++ )
        if ( If_CluHasVar( t, nVarsAll, v ) )
        {
            uSupp |= (1 << v);
            if ( pSupp )
                pSupp[iVar] = pSupp[v];
            iVar++;
        }
    if ( pnVars )
        *pnVars = iVar;
    if ( If_CluSuppIsMinBase( uSupp ) )
        return 0;
    If_CluTruthShrink( t, iVar, nVarsAll, uSupp );
    return 1;
}

// returns the best group found
If_Grp_t If_CluCheck( void * p, int nLutSize, word * pTruth0, int nVars, int iVarStart, int iVarStop, int nLutLeaf, int nLutRoot, 
                     If_Grp_t * pR, word * pFunc0, word * pFunc1, word * pLeftOver, int fHashing )
{
//    int fEnableHashing = 0;
    If_Grp_t G1 = {0}, R = {0};
    unsigned * pHashed = NULL;
    word Truth, pTruth[CLU_WRD_MAX], pF[CLU_WRD_MAX];//, pG[CLU_WRD_MAX];
    int V2P[CLU_VAR_MAX+2], P2V[CLU_VAR_MAX+2], pCanonPerm[CLU_VAR_MAX];
    int i, nSupp, uCanonPhase;
    //int nLutSize = p ? p->pPars->nLutSize : nVars;
    assert( nVars <= CLU_VAR_MAX );
    assert( nVars <= nLutLeaf + nLutRoot - 1 );

    if ( pR )
    {
        pR->nVars = 0;
        *pFunc0 = 0;
        *pFunc1 = 0;
    }

    // canonicize truth table
    //If_CluCopy( pTruth, pTruth0, nLutSize );
    If_CluCopy( pTruth, pTruth0, nVars );    


    // check minnimum base
    If_CluCopy( pF, pTruth, nVars );
    for ( i = 0; i < nVars; i++ )
        V2P[i] = P2V[i] = i;
    // check support
    nSupp = If_CluSupport( pF, nVars );
//Extra_PrintBinary( stdout, &nSupp, 16 );  printf( "\n" );
    if ( !nSupp || !If_CluSuppIsMinBase(nSupp) )
    {
//        assert( 0 );     
        return G1;
    }

    // update the variable order so that the first var was the last one
    if ( iVarStop )
        If_CluMoveVar( pF, nVars, V2P, P2V, 0, nVars-1 );

    if ( G1.nVars == 0 ) 
    {
        s_Count2++;

        // detect easy cofs
        if ( iVarStart == 0 )
            G1 = If_CluDecUsingCofs( pTruth, nVars, nLutLeaf );
        if ( G1.nVars == 0 )
        {
            // perform testing
            G1 = If_CluFindGroup( pF, nVars, iVarStart, iVarStop, V2P, P2V, nLutLeaf, nLutLeaf + nLutRoot == nVars + 1 );
    //        If_CluCheckPerm( pTruth, pF, nVars, V2P, P2V );
            if ( G1.nVars == 0 )
            {
                // perform testing with a smaller set
                if ( nVars < nLutLeaf + nLutRoot - 2 )
                {
                    nLutLeaf--;
                    G1 = If_CluFindGroup( pF, nVars, iVarStart, iVarStop, V2P, P2V, nLutLeaf, nLutLeaf + nLutRoot == nVars + 1 );
                    nLutLeaf++;
                }
                // perform testing with a smaller set
                if ( nLutLeaf > 4 && nVars < nLutLeaf + nLutRoot - 3 )
                {
                    nLutLeaf--;
                    nLutLeaf--;
                    G1 = If_CluFindGroup( pF, nVars, iVarStart, iVarStop, V2P, P2V, nLutLeaf, nLutLeaf + nLutRoot == nVars + 1 );
                    nLutLeaf++;
                    nLutLeaf++;
                }
                if ( G1.nVars == 0 )
                {
                    // perform testing with a different order
                    If_CluReverseOrder( pF, nVars, V2P, P2V, iVarStart );
                    G1 = If_CluFindGroup( pF, nVars, iVarStart, iVarStop, V2P, P2V, nLutLeaf, nLutLeaf + nLutRoot == nVars + 1 );

                    // check permutation
    //                If_CluCheckPerm( pTruth, pF, nVars, V2P, P2V );
                    if ( G1.nVars == 0 )
                    {
                        // remember free set, just in case
//                        for ( i = 0; i < nVars - nLutLeaf; i++ )
///                           G1.pVars[nLutLeaf+i] = P2V[i];
                        // if <XY>, this will not be used
                        // if <XYZ>, this will not be hashed

    /*
                        if ( nVars == 6 )
                        {
                            Extra_PrintHex( stdout, (unsigned *)pF, nVars );  printf( "    " );
                            Kit_DsdPrintFromTruth( (unsigned*)pF, nVars );  printf( "\n" );
                            if ( !If_CutPerformCheck07( (unsigned *)pF, 6, 6, NULL ) )
                                printf( "no\n" );
                        } 
    */
                        if ( pHashed )
                            *pHashed = If_CluGrp2Uns( &G1 );
                        return G1;
                    }
                }
            }
        }
    }

    // derive
    if ( pR )
    {
        int iNewPos;

        If_CluMoveGroupToMsb( pF, nVars, V2P, P2V, &G1 );
        if ( G1.nMyu == 2 )
        {
            Truth = If_CluDeriveDisjoint( pF, nVars, V2P, P2V, &G1, &R );
            iNewPos = R.nVars - 1;
        }
        else
        {
            Truth = If_CluDeriveNonDisjoint( pF, nVars, V2P, P2V, &G1, &R );
            iNewPos = R.nVars - 2;
        }
        assert( R.pVars[iNewPos] == nVars );

        // adjust the functions
        Truth = If_CluAdjust( Truth, G1.nVars );
        if ( R.nVars < 6 )
            pF[0] = If_CluAdjust( pF[0], R.nVars );

//        Kit_DsdPrintFromTruth( (unsigned*)&Truth, G1.nVars ); printf( "  ...1\n" );
//        Kit_DsdPrintFromTruth( (unsigned*)pF, R.nVars );      printf( "  ...1\n" );

        // update the variable order of R so that the new var was the first one
//        if ( iVarStart == 0 )
        {
            int k, V2P2[CLU_VAR_MAX+2], P2V2[CLU_VAR_MAX+2];
            assert( iNewPos >= iVarStart );
            for ( k = 0; k < R.nVars; k++ )
                V2P2[k] = P2V2[k] = k;
            If_CluMoveVar( pF, R.nVars, V2P2, P2V2, iNewPos, iVarStart );
            for ( k = iNewPos; k > iVarStart; k-- )
                R.pVars[k] = R.pVars[k-1];
            R.pVars[iVarStart] = nVars;
        }

//        Kit_DsdPrintFromTruth( (unsigned*)pF, R.nVars ); printf( "  ...2\n" );

        if ( pLeftOver )
        {
            if ( R.nVars < 6 )
                *pLeftOver = If_CluAdjust( pF[0], R.nVars );
            else
                If_CluCopy( pLeftOver, pF, R.nVars );
            If_CluAdjustBig( pLeftOver, R.nVars, nLutSize );
        }

        // perform checking
        if ( 0 )
        {
            If_CluCheckGroup( pTruth, nVars, &G1 );
            If_CluVerify( pTruth, nVars, &G1, &R, Truth, pF );
        } 

        // save functions
        *pR = R;
        if ( pFunc0 )
            *pFunc0 = pF[0];
        if ( pFunc1 )
            *pFunc1 = Truth;
    }

    if ( pHashed )
        *pHashed = If_CluGrp2Uns( &G1 );
    return G1;
}

// returns the best group found
If_Grp_t If_CluCheck3( void * p, int nLutSize, word * pTruth0, int nVars, int nLutLeaf, int nLutLeaf2, int nLutRoot, 
                      If_Grp_t * pR, If_Grp_t * pG2, word * pFunc0, word * pFunc1, word * pFunc2 )
{
//    static int Counter = 0;
    word pLeftOver[CLU_WRD_MAX], Func0, Func1, Func2;
    If_Grp_t G1 = {0}, G2 = {0}, R = {0}, R2 = {0};
    int i;
//    Counter++;
//    if ( Counter == 37590 )
//    {
//        int ns = 0;
//    }

    // check two-node decomposition
    G1 = If_CluCheck( p, nLutSize, pTruth0, nVars, 0, 0, nLutLeaf, nLutRoot + nLutLeaf2 - 1, &R2, &Func0, &Func1, pLeftOver, 0 );
    // decomposition does not exist
    if ( G1.nVars == 0 )
    {
        // check for decomposition with two outputs
        if ( (G1.nMyu == 3 || G1.nMyu == 4) && nLutLeaf == nLutLeaf2 && nVars - nLutLeaf + 2 <= nLutRoot )
        {
            int V2P[CLU_VAR_MAX+2], P2V[CLU_VAR_MAX+2];
            word Func0, Func1, Func2;
            int iVar0, iVar1;

            G1.nVars = nLutLeaf;
            If_CluCopy( pLeftOver, pTruth0, nVars );
            for ( i = 0; i < nVars; i++ )
                V2P[i] = P2V[i] = i;

            If_CluMoveGroupToMsb( pLeftOver, nVars, V2P, P2V, &G1 );
            If_CluDeriveDisjoint4( pLeftOver, nVars, V2P, P2V, &G1, &R, &Func1, &Func2 );

            // move the two vars to the front
            for ( i = 0; i < R.nVars; i++ )
                V2P[i] = P2V[i] = i;
            If_CluMoveVar( pLeftOver, R.nVars, V2P, P2V, R.nVars-2, 0 );
            If_CluMoveVar( pLeftOver, R.nVars, V2P, P2V, R.nVars-1, 1 );
            iVar0 = R.pVars[R.nVars-2];
            iVar1 = R.pVars[R.nVars-1];
            for ( i = R.nVars-1; i > 1; i-- )
                R.pVars[i] = R.pVars[i-2];
            R.pVars[0] = iVar0;
            R.pVars[1] = iVar1;

            Func0 = pLeftOver[0];
            If_CluVerify3( pTruth0, nVars, &G1, &G1, &R, Func1, Func2, Func0 );
            if ( pFunc1 && pFunc2 )
            {
                *pFunc0 = Func0;
                *pFunc1 = Func1;
                *pFunc2 = Func2;
                *pG2 = G1;
                *pR = R;
            }
//                Kit_DsdPrintFromTruth( (unsigned*)pTruth0, nVars );  printf( "\n" );
//                If_CluPrintGroup( &G1 );
            return G1;
        }
        return G1;
    }
    // decomposition exists and sufficient
    if ( R2.nVars <= nLutRoot )
    {
        if ( pG2 )     *pG2 = G2;
        if ( pR )      *pR  = R2;
        if ( pFunc0 )  *pFunc0 = Func0;
        if ( pFunc1 )  *pFunc1 = Func1;
        if ( pFunc2 )  *pFunc2 = 0;
        return G1;
    }

    int nStructType = 1; // LUT structure type
    // the new variable is at the bottom - skip it (iVarStart = 1)
    if ( nStructType == 0 ) // allowed
        G2 = If_CluCheck( p, nLutSize, pLeftOver, R2.nVars, 0, 0, nLutLeaf2, nLutRoot, &R, &Func0, &Func2, NULL, 0 );
    else if ( nStructType == 1 ) // not allowed
        G2 = If_CluCheck( p, nLutSize, pLeftOver, R2.nVars, 1, 0, nLutLeaf2, nLutRoot, &R, &Func0, &Func2, NULL, 0 );
    else if ( nStructType == 2 ) // required
        G2 = If_CluCheck( p, nLutSize, pLeftOver, R2.nVars, 0, 1, nLutLeaf2, nLutRoot, &R, &Func0, &Func2, NULL, 0 );
    else assert( 0 );
    if ( G2.nVars == 0 )
        return G2;

    // remap variables
    for ( i = 0; i < G2.nVars; i++ )
    {
        assert( G2.pVars[i] < R2.nVars );
        G2.pVars[i] = R2.pVars[ (int)G2.pVars[i] ];
    }
    // remap variables
    for ( i = 0; i < R.nVars; i++ )
    {
        if ( R.pVars[i] == R2.nVars )
            R.pVars[i] = nVars + 1;
        else
            R.pVars[i] = R2.pVars[ (int)R.pVars[i] ];
    }

    // decomposition exist
    if ( pG2 )     *pG2 = G2;
    if ( pR )      *pR  = R;
    if ( pFunc0 )  *pFunc0 = Func0;
    if ( pFunc1 )  *pFunc1 = Func1;
    if ( pFunc2 )  *pFunc2 = Func2;

    // verify
    If_CluVerify3( pTruth0, nVars, &G1, &G2, &R, Func1, Func2, Func0 );
    return G1;
}

////////////////////////////////////////////////////
//   Support minimization added on Sep 8, 2023    //
////////////////////////////////////////////////////

static word s_Truths6[6] = {
    ABC_CONST(0xAAAAAAAAAAAAAAAA),
    ABC_CONST(0xCCCCCCCCCCCCCCCC),
    ABC_CONST(0xF0F0F0F0F0F0F0F0),
    ABC_CONST(0xFF00FF00FF00FF00),
    ABC_CONST(0xFFFF0000FFFF0000),
    ABC_CONST(0xFFFFFFFF00000000)
};

static word s_Truths6Neg[6] = {
    ABC_CONST(0x5555555555555555),
    ABC_CONST(0x3333333333333333),
    ABC_CONST(0x0F0F0F0F0F0F0F0F),
    ABC_CONST(0x00FF00FF00FF00FF),
    ABC_CONST(0x0000FFFF0000FFFF),
    ABC_CONST(0x00000000FFFFFFFF)
};

static word s_PMasks[5][3] = {
    { ABC_CONST(0x9999999999999999), ABC_CONST(0x2222222222222222), ABC_CONST(0x4444444444444444) },
    { ABC_CONST(0xC3C3C3C3C3C3C3C3), ABC_CONST(0x0C0C0C0C0C0C0C0C), ABC_CONST(0x3030303030303030) },
    { ABC_CONST(0xF00FF00FF00FF00F), ABC_CONST(0x00F000F000F000F0), ABC_CONST(0x0F000F000F000F00) },
    { ABC_CONST(0xFF0000FFFF0000FF), ABC_CONST(0x0000FF000000FF00), ABC_CONST(0x00FF000000FF0000) },
    { ABC_CONST(0xFFFF00000000FFFF), ABC_CONST(0x00000000FFFF0000), ABC_CONST(0x0000FFFF00000000) }
};

static word s_PPMasks[5][6][3] = {
    { 
        { ABC_CONST(0x0000000000000000), ABC_CONST(0x0000000000000000), ABC_CONST(0x0000000000000000) }, // 0 0  
        { ABC_CONST(0x9999999999999999), ABC_CONST(0x2222222222222222), ABC_CONST(0x4444444444444444) }, // 0 1  
        { ABC_CONST(0xA5A5A5A5A5A5A5A5), ABC_CONST(0x0A0A0A0A0A0A0A0A), ABC_CONST(0x5050505050505050) }, // 0 2 
        { ABC_CONST(0xAA55AA55AA55AA55), ABC_CONST(0x00AA00AA00AA00AA), ABC_CONST(0x5500550055005500) }, // 0 3 
        { ABC_CONST(0xAAAA5555AAAA5555), ABC_CONST(0x0000AAAA0000AAAA), ABC_CONST(0x5555000055550000) }, // 0 4 
        { ABC_CONST(0xAAAAAAAA55555555), ABC_CONST(0x00000000AAAAAAAA), ABC_CONST(0x5555555500000000) }  // 0 5 
    },
    { 
        { ABC_CONST(0x0000000000000000), ABC_CONST(0x0000000000000000), ABC_CONST(0x0000000000000000) }, // 1 0  
        { ABC_CONST(0x0000000000000000), ABC_CONST(0x0000000000000000), ABC_CONST(0x0000000000000000) }, // 1 1  
        { ABC_CONST(0xC3C3C3C3C3C3C3C3), ABC_CONST(0x0C0C0C0C0C0C0C0C), ABC_CONST(0x3030303030303030) }, // 1 2 
        { ABC_CONST(0xCC33CC33CC33CC33), ABC_CONST(0x00CC00CC00CC00CC), ABC_CONST(0x3300330033003300) }, // 1 3 
        { ABC_CONST(0xCCCC3333CCCC3333), ABC_CONST(0x0000CCCC0000CCCC), ABC_CONST(0x3333000033330000) }, // 1 4 
        { ABC_CONST(0xCCCCCCCC33333333), ABC_CONST(0x00000000CCCCCCCC), ABC_CONST(0x3333333300000000) }  // 1 5 
    },
    { 
        { ABC_CONST(0x0000000000000000), ABC_CONST(0x0000000000000000), ABC_CONST(0x0000000000000000) }, // 2 0  
        { ABC_CONST(0x0000000000000000), ABC_CONST(0x0000000000000000), ABC_CONST(0x0000000000000000) }, // 2 1  
        { ABC_CONST(0x0000000000000000), ABC_CONST(0x0000000000000000), ABC_CONST(0x0000000000000000) }, // 2 2 
        { ABC_CONST(0xF00FF00FF00FF00F), ABC_CONST(0x00F000F000F000F0), ABC_CONST(0x0F000F000F000F00) }, // 2 3 
        { ABC_CONST(0xF0F00F0FF0F00F0F), ABC_CONST(0x0000F0F00000F0F0), ABC_CONST(0x0F0F00000F0F0000) }, // 2 4 
        { ABC_CONST(0xF0F0F0F00F0F0F0F), ABC_CONST(0x00000000F0F0F0F0), ABC_CONST(0x0F0F0F0F00000000) }  // 2 5 
    },
    { 
        { ABC_CONST(0x0000000000000000), ABC_CONST(0x0000000000000000), ABC_CONST(0x0000000000000000) }, // 3 0  
        { ABC_CONST(0x0000000000000000), ABC_CONST(0x0000000000000000), ABC_CONST(0x0000000000000000) }, // 3 1  
        { ABC_CONST(0x0000000000000000), ABC_CONST(0x0000000000000000), ABC_CONST(0x0000000000000000) }, // 3 2 
        { ABC_CONST(0x0000000000000000), ABC_CONST(0x0000000000000000), ABC_CONST(0x0000000000000000) }, // 3 3 
        { ABC_CONST(0xFF0000FFFF0000FF), ABC_CONST(0x0000FF000000FF00), ABC_CONST(0x00FF000000FF0000) }, // 3 4 
        { ABC_CONST(0xFF00FF0000FF00FF), ABC_CONST(0x00000000FF00FF00), ABC_CONST(0x00FF00FF00000000) }  // 3 5 
    },
    { 
        { ABC_CONST(0x0000000000000000), ABC_CONST(0x0000000000000000), ABC_CONST(0x0000000000000000) }, // 4 0  
        { ABC_CONST(0x0000000000000000), ABC_CONST(0x0000000000000000), ABC_CONST(0x0000000000000000) }, // 4 1  
        { ABC_CONST(0x0000000000000000), ABC_CONST(0x0000000000000000), ABC_CONST(0x0000000000000000) }, // 4 2 
        { ABC_CONST(0x0000000000000000), ABC_CONST(0x0000000000000000), ABC_CONST(0x0000000000000000) }, // 4 3 
        { ABC_CONST(0x0000000000000000), ABC_CONST(0x0000000000000000), ABC_CONST(0x0000000000000000) }, // 4 4 
        { ABC_CONST(0xFFFF00000000FFFF), ABC_CONST(0x00000000FFFF0000), ABC_CONST(0x0000FFFF00000000) }  // 4 5 
    }    
};

#define ABC_SWAP(Type, a, b)  { Type t = a; a = b; b = t; }

static inline word Abc_Tt6SwapVars( word t, int iVar, int jVar )
{
    word * s_PMasks = s_PPMasks[iVar][jVar];
    int shift = (1 << jVar) - (1 << iVar);
    assert( iVar < jVar );
    return (t & s_PMasks[0]) | ((t & s_PMasks[1]) << shift) | ((t & s_PMasks[2]) >> shift);
}
static inline void Abc_TtSwapVars( word * pTruth, int nVars, int iVar, int jVar )
{
    if ( iVar == jVar )
        return;
    if ( jVar < iVar )
        ABC_SWAP( int, iVar, jVar );
    assert( iVar < jVar && jVar < nVars );
    if ( nVars <= 6 )
    {
        pTruth[0] = Abc_Tt6SwapVars( pTruth[0], iVar, jVar );
        return;
    }
    if ( jVar <= 5 )
    {
        word * s_PMasks = s_PPMasks[iVar][jVar];
        int nWords = Abc_TtWordNum(nVars);
        int w, shift = (1 << jVar) - (1 << iVar);
        for ( w = 0; w < nWords; w++ )
            pTruth[w] = (pTruth[w] & s_PMasks[0]) | ((pTruth[w] & s_PMasks[1]) << shift) | ((pTruth[w] & s_PMasks[2]) >> shift);
        return;
    }
    if ( iVar <= 5 && jVar > 5 )
    {
        word low2High, high2Low;
        word * pLimit = pTruth + Abc_TtWordNum(nVars);
        int j, jStep = Abc_TtWordNum(jVar);
        int shift = 1 << iVar;
        for ( ; pTruth < pLimit; pTruth += 2*jStep )
            for ( j = 0; j < jStep; j++ )
            {
                low2High = (pTruth[j] & s_Truths6[iVar]) >> shift;
                high2Low = (pTruth[j+jStep] << shift) & s_Truths6[iVar];
                pTruth[j] = (pTruth[j] & ~s_Truths6[iVar]) | high2Low;
                pTruth[j+jStep] = (pTruth[j+jStep] & s_Truths6[iVar]) | low2High;
            }
        return;
    }
    {
        word * pLimit = pTruth + Abc_TtWordNum(nVars);
        int i, iStep = Abc_TtWordNum(iVar);
        int j, jStep = Abc_TtWordNum(jVar);
        for ( ; pTruth < pLimit; pTruth += 2*jStep )
            for ( i = 0; i < jStep; i += 2*iStep )
                for ( j = 0; j < iStep; j++ )
                    ABC_SWAP( word, pTruth[iStep + i + j], pTruth[jStep + i + j] );
        return;
    }    
}
static inline int Abc_TtMinBase( word * pTruth, int * pVars, int nVars, int nVarsAll ) 
{
    int i, k;
    assert( nVars <= nVarsAll );
    for ( i = k = 0; i < nVars; i++ )
    {
        if ( !Abc_TtHasVar( pTruth, nVarsAll, i ) )
            continue;
        if ( k < i )
        {
            if ( pVars ) pVars[k] = pVars[i];
            Abc_TtSwapVars( pTruth, nVarsAll, k, i );
        }
        k++;
    }
    if ( k == nVars )
        return k;
    assert( k < nVars );
//    assert( k == Abc_TtSupportSize(pTruth, nVars) );
    return k;
}

If_Grp_t If_CluCheckTest( int Size, int nLutSize, word * Truth, int nLeaves, 
                          If_Grp_t * pR, If_Grp_t * pG2, word * pFunc0, word * pFunc1, word * pFunc2, int * pnVarsNew, int * pVarPerm )
{
        char * pStr = Size == 2 ? (char *)"66" : (Size == 3 ? (char *)"666" : NULL);
        int i, nLutLeaf2, nLutRoot;
        int nLutLeaf, nRootLeaf;
        If_Grp_t G1 = {0};
        assert( nLeaves > 6 );

        ////////////////////////////////////////////////////////////
        // perform support minimization
        int nLeavesOld = nLeaves;
        for ( i = 0; i < nLeaves; i++ )
            pVarPerm[i] = i;           
        nLeaves = Abc_TtMinBase( Truth, pVarPerm, nLeaves, nLeaves );
        if ( nLeaves < nLeavesOld ) {
            printf( "Support minimization reduced %d variables. New support is ", nLeavesOld - nLeaves );
            printf( "{" );
            for ( i = 0; i < nLeaves; i++ )
                printf( " %d", pVarPerm[i] );
            printf( " }\n" );
        }
        // return the number of variables after minimization
        *pnVarsNew = nLeaves;
        if ( nLeaves <= 6 ) {
            printf( "The support does not exceed the LUT size. Decomposition is not performed.\n" );
            return G1;
        }
        // after support minimization, the number of input variables may be reduced from nLeavesOld to nLeaves
        // while the indexes of the remaining nLeaves variables are listed in the array pVarPerm[]
        ////////////////////////////////////////////////////////////

        // make sure the function is support-minimal
        for ( int v = 0; v < nLeaves; v++ )
            if ( !Abc_TtHasVar( Truth, nLeaves, v ) )
            {
                printf( "The function is non-support minimal. Decomposition is not performed.\n" );
                return G1;
            }

        // quit if parameters are wrong
        int Length = strlen(pStr);
        if ( Length != 2 && Length != 3 )
        {
            printf( "Wrong LUT struct (%s)\n", pStr );
            return G1;
        }
        for ( i = 0; i < Length; i++ )
            if ( pStr[i] - '0' < 3 || pStr[i] - '0' > 6 )
            {
                printf( "The LUT size (%d) should belong to {3,4,5,6}.\n", pStr[i] - '0' );
                return G1;
            }

        nLutLeaf  =                   pStr[0] - '0';
        nLutLeaf2 = ( Length == 3 ) ? pStr[1] - '0' : 0;
        nLutRoot  =                   pStr[Length-1] - '0';
        if ( nLeaves > nLutLeaf - 1 + (nLutLeaf2 ? nLutLeaf2 - 1 : 0) + nLutRoot )
        {
            printf( "The cut size (%d) is too large for the LUT structure %s.\n", nLeaves, pStr );
            return G1;
        }
        // consider easy case
        if ( nLeaves <= Abc_MaxInt( nLutLeaf2, Abc_MaxInt(nLutLeaf, nLutRoot) ) )
            return G1;

        // derive the first group
        if ( Length == 2 )
            G1 = If_CluCheck( NULL, nLutSize, Truth, nLeaves, 0, 0, nLutLeaf, nLutRoot, pR, pFunc0, pFunc1, NULL, 0 );
        else
            G1 = If_CluCheck3( NULL, nLutSize, Truth, nLeaves, nLutLeaf, nLutLeaf2, nLutRoot, pR, pG2, pFunc0, pFunc1, pFunc2 );
        return G1;
}

/*************************************************************
                  Reading input data
**************************************************************/

static inline int Abc_Base2Log( unsigned n )   { int r; if ( n < 2 ) return (int)n; for ( r = 0, n--; n; n >>= 1, r++ ) {}; return r; }
static inline int Abc_Hex2Int( char c ) {
    int Digit = 0;
    if ( c >= '0' && c <= '9' )
        Digit = c - '0';
    else if ( c >= 'A' && c <= 'F' )
        Digit = c - 'A' + 10;
    else if ( c >= 'a' && c <= 'f' )
        Digit = c - 'a' + 10;
    else assert( 0 );
    assert( Digit >= 0 && Digit < 16 );
    return Digit;
}
static inline void Abc_PrintHexTruth( word * Truth, int nVars )  {
    int k, Digit, nDigits = Abc_TtHexDigitNum(nVars);
    for ( k = nDigits-1; k >= 0; k-- ) {
        Digit = (int)((Truth[k/16] >> ((k%16) * 4)) & 15);
        if ( Digit < 10 )
            printf( "%d", Digit );
        else
            printf( "%c", 'A' + Digit-10 );
    }
}
void Kit_DsdPrintFromTruth( unsigned * pTruth, int nVars )
{
    Abc_PrintHexTruth( (word *)pTruth, nVars );
}
static inline int Abc_ReadHexTruth( char * pInput, word * Truth )
{
    int nChars = strlen(pInput), i; 
    int nVars  = Abc_Base2Log(4*nChars);
    if ( (1 << nVars) != 4*nChars ) {
        printf( "The input string length (%d chars) does not match the size (%d bits) of the truth table of %d-var function.\n", 
            nChars, 1<<nVars, nVars );
        return 0;
    }
    word Num = 0;
    for ( i = nChars-1; i >= 0; i-- ) {
        Num |= (word)Abc_Hex2Int(pInput[nChars-1-i]) << ((i & 0xF) * 4);
        if ( (i & 0xF) == 0 )
            Truth[i>>4] = Num, Num = 0;
    }
    assert( Num == 0 );
    if ( nVars < 6 ) 
        Truth[0] = Abc_Tt6Stretch( Truth[0], nVars );
    if ( nVars < 7 )
        Truth[1] = Truth[0];
    printf( "Finished entring %d-input function: ", nVars );
    Abc_PrintHexTruth( Truth, nVars );
    printf( "\n" );
    return nVars;
}

/*************************************************************
                    main() procedure
**************************************************************/

// int main(int argc, char ** argv)
// {
//     if ( argc == 1 )
//     {
//         printf( "usage:  %s <truth_table>\n", argv[0] );
//         printf( "                  this program decomposes an N-input function into two or three 6-LUTs (S66 or S666)\n" );
//         printf( "  <truth_table> : the truth table of the function in the hexadecimal notation\n" );
//         return 1;
//     }
//     else
//     {
//         int nLutSize = 6;
//         word Truth[CLU_WRD_MAX] = {0};
//         int nVars = Abc_ReadHexTruth( argv[1], Truth );
//         if ( nVars <= 6 ) {
//             printf( "Function has less than 6 inputs.\n");
//             return 1;
//         }

//         word Func0, Func1, Func2;
//         If_Grp_t G1 = {0}, G2 = {0}, R = {0};
//         int nVarsNew = nVars;            // the number of variables afer support minimization
//         int pVarPerm[CLU_VAR_MAX] = {0}; // the remaining variables after support minimization
//         G1 = If_CluCheckTest( 2, nLutSize, Truth, nVars, &R, &G2, &Func0, &Func1, &Func2, &nVarsNew, pVarPerm );
//         if ( G1.nVars > 0 ) {
//             If_CluPrintConfig( nVarsNew, &G1, &R, Func1, &Func0 );
//             //If_CluPrintGroup( &G1 );
//             //Kit_DsdPrintFromTruth( (unsigned *)&Func1, 6 ); printf( "\n" );
//             //Kit_DsdPrintFromTruth( (unsigned *)&Func0, 6 ); printf( "\n" );
//         }
//         else
//             printf( "This function cannot be decomposed into S66.\n" );

//         return 1;
//     }
// }

/*************************************************************
                     end of file
**************************************************************/
