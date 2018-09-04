/*
 * Test rnamoves.c
 * Alternative move set definitions for RNAs
 * Use the Check Unit Testing Framework
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>  /* to implement assign() */

#include <check.h>

#include "rnamoves.h"
#include "rnamoves.c"   /* to test non-exported functions. Yea, it's dirty. */
#include "stapel.h"


/*************************
 **  Utility Functions  **
 *************************/

/**
 * @brief Check two arrays of NULL-terminated char strings for deep equality
 *          using Check assertions. Check second array in reversed order.
 * @param First string array
 * @param Size of first string array
 * @param Second string array
 * @param Size of second string array
 */
void deep_eq_str_ary( char** a, int size_a, char** b, int size_b);

/**
 * @brief Check two arrays of integers for deep equality
 *          using Check assertions.
 * @param First int array
 * @param Size of first int array
 * @param Second int array
 * @param Size of second int array
 */
void deep_eq_int_ary( int* a, int size_a, int* b, int size_b);

/**
 * @brief Pop all structures off the global structure stack into an array
 * @param size Pointer to an integer in which the number of structures found
 *          is stored
 * @return Array of strings containing all structures found on the global
 *          structure stack
 */
char** pop_stapel( int* size);

/**
 * @brief Assign the passed values to the passed array.
 * @param Integer array
 * @param nparam Number of values to assign.
 * @param ... Values to assign to the array
 */
void assign_int( int* array, int nparam, ...);

/**
 * @brief Assign the passed values to the passed array.
 * @param Array for strings, i.e. zero-terminated char arrays
 * @param nparam Number of values to assign.
 * @param ... Values to assign to the array
 */
void assign_str( char** array, int nparam, ...);

/**
 * @brief Test RNA2_move_it by comparing the structures it generates with a
 *          passed list of structures that it SHOULD generate
 * @param seq Sequence
 * @param str Structure as dot-bracket string
 * @param noLP flag,set 1 to avoid lonely pairs
 * @param shift flag,set 1 to allow shift moves
 * @param nb_size actual number of neighbor structures
 * @param ... neighbor structures in the order they are created by RNA2_move_it
 */
void test_RNA2_move_it( char* seq, char* str, int shift, int noLP,
                        int nb_size, ...);

/**
 * @brief A simple in-situ QuickSort implementation for C strings. Use
 *  sort_str() for a more convenient interface.
 * @param ary Array of C-style strings, i.e. zero-terminated char arrays
 * @param l Left index of range to sort
 * @param r Right index of range to sort
 */

void quickSort( char* a[], int l, int r);
/**
 * @brief QuickSort an array of zero-terminated strings
 * @param ary Array of strings (i.e. array of char*)
 * @param size Number of strings in array
 */
inline void sort_str( char* ary[], int size)
{
    quickSort( ary, 0, size-1);
}


/*********************************
 **  Start of test definitions  **
 *********************************/

/* Test sort_str() */
START_TEST( test_sort_str )
{
    int i;
    char* strings[] =
    {
        "Foo",
        "Baz",
        "Zap",
        "Bar",
        "Aww",
        "aww"
    };
    int size = 6;
    char* sorted[] =
    {
        "Aww",
        "Bar",
        "Baz",
        "Foo",
        "Zap",
        "aww"
    };
    sort_str( strings, size);
    deep_eq_str_ary(strings, size, sorted, size);
}
END_TEST

/* Check for correct value of external MYTURN variable */
START_TEST ( test_myturn_val )
{
    ck_assert_int_eq( MYTURN, 4);
}
END_TEST

/* Test the rnamoves_struct data structure */
START_TEST ( test_rnamoves_struct)
{
    char* seq = "AAAUUUGGGCCC";
    char* str = "............";
    int shift = 1, noLP = 0;
    RNA2_init( seq, shift, noLP);
    ck_assert_str_eq( RNA2_get_seq(), seq);
    ck_assert_str_eq( RNA2_get_form(), str);
    ck_assert_int_eq( r2d->shift, shift);
    ck_assert_int_eq( r2d->noLP, noLP);
    /*ck_assert_int_eq( 1, 0);    /* THOU SHALL FAIL! */

    RNA2_free();
}
END_TEST

/* Test generation of pair table from dot-bracket strings */
START_TEST( test_parse_structure )
{
    char* seq =    "tTuAaAGgGCcC";
    char* str =    "............";
    int str_pt[] = { UNPRD, UNPRD, UNPRD,
                      UNPRD, UNPRD, UNPRD,
                      UNPRD, UNPRD, UNPRD,
                      UNPRD, UNPRD, UNPRD
                    };
    RNA2_init( seq, 1, 1);
    strcpy(r2d->form, str);
    parse_structure();
    deep_eq_int_ary( r2d->pt, r2d->len, str_pt, r2d->len);

    /*        012345678901  */
    str =    "((...))(...)";
    assign_int( str_pt, 12,
                6, 5, UNPRD,
                UNPRD, UNPRD, 1,
                0, 11, UNPRD,
                UNPRD, UNPRD, 7
                );
    strcpy(r2d->form, str);
    parse_structure();
    deep_eq_int_ary( r2d->pt, r2d->len, str_pt, r2d->len);

    RNA2_free();
}
END_TEST

/* Test loneliness tests like isLonely(), is_ins_lonely(), insGrowLonely() etc */
START_TEST( test_loneliness_tests)
{
    char* seq = "GGGGGUGCCCCC";
    char* str[] = {
        ".(((....))).", "..((....))..", "..(......)..", ".((......).)",
        ".((.(...)).)", ".(((..()))).", "(((......)))", "((((.)(.))))"
    };
    int str_size = 8;
    int insIsLone[] = { 0, 0, 1, 1, 1, 0, 1, 1};
    int outIsLone[] = { 0, 1, 1, 1, 1, 0, 0, 0};
    int isLone[] = { 0, 0, 1, 1, 1, 0, 0, 0};
    int insGrowLone[] = { 1, 1, 0, 0, 0, 1, 0, 0};
    int outGrowLone[] = { 1, 0, 0, 0, 0, 1, 0, 0};
    int i;
    char* oldForm;
    RNA2_init( seq, 0, 0);
    oldForm = r2d->form;
    for( i=0; i<str_size; i++)
    {
        r2d->form = str[i];
        parse_structure();
        ck_assert_int_eq( is_ins_lonely(2, 9), insIsLone[i]);
        ck_assert_int_eq( is_out_lonely(2, 9), outIsLone[i]);
        ck_assert_int_eq( isLonely(2, 9), isLone[i]);
        ck_assert_int_eq( insGrowLonely(2, 9), insGrowLone[i]);
        ck_assert_int_eq( outGrowLonely(2, 9), outGrowLone[i]);
    }
    r2d->form = oldForm;
    RNA2_free();
}
END_TEST

/* Test generation of neighbor structures, without shift moves, with lonely pairs */
START_TEST( test_rnamoves_noshift_lp)
{
    int noLP = 0;
    int shift = 0;

    char* seq = "GGGGUUUCCCC";
    char* str = "((((...))))";
    int nb_size = 4;
    test_RNA2_move_it( seq, str, shift, noLP, nb_size,
                       "(((.....)))",
                       "((.(...).))",
                       "(.((...)).)",
                       ".(((...)))."
                       );

    seq = "GGUUUCC";
    str = ".......";
    nb_size = 5;                /* Number of passed neighbors */
    test_RNA2_move_it( seq, str, shift, noLP, nb_size,
                       ".(....)",
                       ".(...).",
                       "(.....)",
                       "(....).",
                       "(...).."
                       );
}
END_TEST

/* Test generation of neighbor structures, with shift moves, with lonely pairs */
START_TEST( test_rnamoves_shift_lp)
{
    int shift = 1;
    int noLP = 0;

    char* seq = "GGGGGGGGGGGGAAAUUUUUUUGUUUU";
    char* str = ".(((....(((.....)))....))).";

    int nb_size = 40;
    test_RNA2_move_it( seq, str, shift, noLP, nb_size,
				".(((...((((.....))))...))).", ".(((.(..(((.....)))..).))).",
				"((((....(((.....)))....))))", ".(((....((((...))))....))).",
				"..((....(((.....)))....))..", ".(((.....((.....)).....))).",
				".(((....((.......))....))).", ".(((..(.(((.....))).)..))).",
				".((.....(((.....))).....)).", ".(.(....(((.....)))....).).",
				".(((....(.(.....).)....))).", ".(((....(((....).))....))).",
				".(((.(..(((.....))).)..))).", ".(((..(.(((.....)))..).))).",
				".(((....((.(....)))....))).", ".(((..(.(((.....))))...))).",
				".(((...((((.....))).)..))).", ".(((....(((.....)))....)).)",
				".((.(...(((.....)))....))).", ".(((...(.((.....)))....))).",
				".(((....(((.....)).)...))).", "(.((....(((.....)))....))).",
				".((((...(((.....)))..).))).", ".((((...(((.....))).)..))).",
				".(((.(..(((.....))))...))).", ".(((...((((.....)))..).))).",
				".(((....((..(...)))....))).", ".((..(..(((.....)))....))).",
				".(((..(..((.....)))....))).", ".(((....(((.....))..)..))).",
				".((((...(((.....))))...))).", ".(((....(((.....)))..)..)).",
				".((....((((.....)))....))).", ".(((....(((.....))))....)).",
				".((((....((.....)))....))).", ".((...(.(((.....)))....))).",
				".(((....(((.....))).)...)).", ".(((.(...((.....)))....))).",
				".(((....(((.....))...).))).", ".(((.....((.....))(...))))."
                       );
}
END_TEST

/* Test generation of neighbor structures, without shift moves, without lonely pairs */
START_TEST( test_rnamoves_noshift_nolp )
{
    int shift = 0;
    int noLP = 1;

    char* seq = "GGGGUUUCCCC";
    char* str = "((((...))))";

    int nb_size = 2;
    test_RNA2_move_it( seq, str, shift, noLP, nb_size,
                       "(((.....)))",
                       ".(((...)))."
                       );

    str = "...........";
    nb_size = 12;
    test_RNA2_move_it( seq, str, shift, noLP, nb_size,
                       "((......)).", ".((....))..", ".((.....)).",
                       "((.....))..", "..((....)).", ".((......))",
                       "..((...))..", "((.......))", "..((.....))",
                       "((....))...", ".((...))...", "((...))...."
                       );
}
END_TEST

/* Test generation of neighbor structures, with shift moves, without lonely pairs */
START_TEST( test_rnamoves_shift_nolp)
{
    int shift = 1;
    int noLP = 1;

    char* seq =        "GGGGGGGGUUUUUUU";
    char* str =        "((..(((...)))))";
//               01234567890123456789012
//               0         1         2
    int nb_size = 4;
    test_RNA2_move_it( seq, str, shift, noLP, nb_size,
                       "....(((...)))..",
                       "((..((.....))))",
                       "(((..((...)))))",
                       "((...((...)).))"
                       );

    seq = "AAGUGAUACCAGCAUCGUCUUGAUGCCCUUGGCAGCACUUCAU";
    str = "(((((((......)))).)))((((((...)))).....))..";
    nb_size = 9;
    test_RNA2_move_it( seq, str, shift, noLP, nb_size,
                       "(((((((......))))).))((((((...)))).....))..",
                       "(((((((......)))).)))(((((.....))).....))..",
                       "(((((((......)))).)))((.(((...)))......))..",
                       "(((((((......)))).)))..((((...)))).........",
                       "((((((........))).)))((((((...)))).....))..",
                       "(((.(((......)))..)))((((((...)))).....))..",
                       "((.((((......))))..))((((((...)))).....))..",
                       ".((((((......)))).))(((((((...)))).....))).",
                       ".((((((......)))).)).((((((...)))).....)).."
                       );
}
END_TEST

/* template */
START_TEST( test_name )
{

}
END_TEST

/* ALSO ADD NEW TESTS TO TEST SUITE IN rnamoves_suite()! */


/***************************
 * End of test definitions *
 ***************************/

/**
 * @brief Test the function RNA2_move_it by comparing the generated neighbor
 *  to a passed list of structures. The passed structures are sorted, so their
 *  order is irrelevant. Results will be reported by the test framework (Check)
 * @param seq Sequence
 * @param str Structure of which to generate neighbors
 * @param shift Flag: If 1, allow shift moves
 * @param noLP Flag: If 1, avoid generation of non-canonical structures (i.e.
 *  structures containing isolated or lonely pairs)
 * @param nb_size Number of passed neighbor structures (cf. '...')
 * @param ... several structure dot-bracket strings to which the generated
 *  neighbors should be compared to. This list will be sorted, so the order
 *  in which structures are passed is irrelevant.
 */
void test_RNA2_move_it( char* seq, char* str, int shift, int noLP, int nb_size, ...)
{
    int i;
    va_list valist; /* Variable args list */
    va_start(valist, nb_size);
    int seqLen = strlen(seq);
    char** gen_nb = NULL;   /* To store generated neighbors */
    int gen_nb_size;        /* Number of generated structures */
    char** nb = (char**) malloc( 3*seqLen*seqLen*sizeof(char*));

    ini_stapel( seqLen);

    for( i=0; i<nb_size; i++)
        nb[i] = va_arg( valist, char*);
    sort_str( nb, nb_size);
    RNA2_init( seq, shift, noLP);

    RNA2_move_it( str);
    gen_nb = pop_stapel( &gen_nb_size);    /* Collect neighbors from stapel */
    sort_str( gen_nb, gen_nb_size);         /* Sort generated neighbors */
    fprintf( stderr, "\"%s\" (input)\n", str);
    for(i = 0; i<gen_nb_size; i++)      /* Output generated neighbors */
         fprintf( stderr, "\"%s\",\n", gen_nb[i]);
    deep_eq_str_ary( gen_nb, gen_nb_size, nb, nb_size);

    free( gen_nb);
    free( nb);
    va_end( valist);
    RNA2_free();
    free_stapel();
}

/**
 * @brief Test two arrays of C-style strings for deep equality. The results
 *  are reported via the testing framework (Check).
 * @param a First array of strings
 * @param size_a Size of first array
 * @param b Second array of strings
 * @param size_b Size of second array
 */
void deep_eq_str_ary( char** a, int size_a, char** b, int size_b)
{
    int i;
    ck_assert_int_eq( size_a, size_b);
    for( i=0; i<size_a; i++)
    {
        /*fprintf( stderr, "%d\n", i); */
        ck_assert_str_eq(a[i], b[i]);
    }

    return;
}

/**
 * @brief Test two arrays of integers for deep equality. The results
 *  are reported via the testing framework (Check).
 * @param a First array of integers
 * @param size_a Size of first array
 * @param b Second array of integers
 * @param size_b Size of second array
 */
void deep_eq_int_ary( int* a, int size_a, int* b, int size_b)
{
    int i;
    ck_assert_int_eq( size_a, size_b);
    for( i=0; i<size_a; i++)
        ck_assert_int_eq(a[i], b[i]);

    return;
}

/**
 * @brief Function that pops all elements from Vienna's global structure stack
 *  and returns them as a single array of strings.
 * @param size Pointer to an integer that is set to the number of popped
 *  structures
 * @return Pointer to an array of structure dot-bracket strings
 */
char** pop_stapel( int* size)
{
    char* cur = NULL;
    int seqLen = RNA2_get_len();
    char** structs = (char**) malloc( 3*seqLen*seqLen*sizeof(char*));
    *size = 0;
    while( cur = pop())
        structs[(*size)++] = cur;

    return structs;
}

/**
 * @brief Utility function to assign a list of integer values to an integer array.
 * @param array The array the values should be assigned to. Needs to be of
 *  sufficient size!
 * @param nparam Number of integers to assign
 * @param ... integers to assign to the array
 */
void assign_int( int* array, int nparam, ...)
{
    va_list valist; /* Variable arguments list */
    int i = 0;

    va_start(valist, nparam);
    for (i = 0; i < nparam; i++)
        array[i] = va_arg(valist, int);

    va_end(valist);
    return;
}

/**
 * @brief Utility function to assign a list of strings to a char* array.
 * @param array The array the values should be assigned to. Needs to be of
 *  sufficient size!
 * @param nparam Number of strings to assign
 * @param ... strings to assign to the array
 */
void assign_str( char** array, int nparam, ...)
{
    va_list valist; /* Variable arguments list */
    int i = 0;

    va_start(valist, nparam);
    for (i = 0; i < nparam; i++)
        array[i] = va_arg(valist, char*);

    va_end(valist);
    return;
}


/************************************************************
 * QuickSort takenfrom
 * http://www.comp.dit.ie/rlawlor/Alg_DS/sorting/quickSort.c
 * Adapted to sort string arrays, core dumps fixed.
 ************************************************************/

int partition( char* a[], int l, int r) {
    int i, j;
    char *pivot, *t;
    pivot = a[l];
    i = l; j = r+1;

    while( 1)
    {
        /* FK: first compare i and r to prevent core dumps */
        do ++i; while( i<=r && strcmp( a[i], pivot) <= 0 );
        do --j; while( strcmp( a[j], pivot) > 0);
        if( i >= j ) break;
        t = a[i]; a[i] = a[j]; a[j] = t;
    }
    t = a[l]; a[l] = a[j]; a[j] = t;
    return j;
}

void quickSort( char* a[], int l, int r)
{
    int j;

    if( l < r )
    {
        // divide and conquer
        j = partition( a, l, r);
        quickSort( a, l, j-1);
        quickSort( a, j+1, r);
    }
}

/***************************************
 * end of QuickSort
 * *************************************/

/***********************
 **  Setup Testsuite  **
 ***********************/

/* Create a test suite
 * Add new tests to suite here!
 * Taken from the Check tutorial.
 */
Suite* rnamoves_suite(void)
{
    Suite *s;
    TCase *tc_core;

    s = suite_create("rnamoves");

    /* Core test case */
    tc_core = tcase_create("Core");

    /* tcase_add_test(tc_core, <test_name> ); */

    tcase_add_test(tc_core, test_sort_str);
    tcase_add_test(tc_core, test_myturn_val);
    tcase_add_test(tc_core, test_rnamoves_struct);
    tcase_add_test(tc_core, test_parse_structure);
    tcase_add_test(tc_core, test_loneliness_tests);
    tcase_add_test(tc_core, test_rnamoves_noshift_lp);
    tcase_add_test(tc_core, test_rnamoves_shift_lp);
    tcase_add_test(tc_core, test_rnamoves_noshift_nolp);
    tcase_add_test(tc_core, test_rnamoves_shift_nolp);

    suite_add_tcase(s, tc_core);
    return s;
}

/* Run test suite */
int main(int argc, char **args)
{
    /* printf( "Hello World! check_rnamoves.c:main()\n"); */
    int number_failed;
    Suite *s;
    SRunner *sr;

    s = rnamoves_suite();
    sr = srunner_create(s);

    srunner_run_all(sr,CK_VERBOSE);/* CK_... SILENT MINIMAL NORMAL VERBOSE etc*/
    number_failed = srunner_ntests_failed(sr);
    srunner_free(sr);
    return (number_failed == 0) ? EXIT_SUCCESS : EXIT_FAILURE;

    return 0;
}
