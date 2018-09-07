#ifndef _rnamoves_h
#define _rnamoves_h
/**
 * Alternative move set definitions for RNAs
 * for use in the RNA landscape analysis tool "barriers"
 *
 * Author(s): Felix Kuehnl, Sebastian Will
 * Date:      Oct-Nov 2015
 *
 * @author Felix Kuehnl, Sebastian Will
 * @date Oct 2015
 */


/**
 * @brief initialize RNA2 move generator
 * @param sequence RNA sequence
 * @param structure of RNA sequence as dot-bracket string
 * @param shift flag, turn on generation of shift moves         
 * @param nolp flag, turn on move generation for noLP structures
 */
void RNA2_init(char const * const sequence, int const shift, int const nolp);


/**
 * @brief generate neighbors according to move set. The flags 'shift'
 *  and 'noLP' of the r2d structure are taken into account.
 * @param struc source structure as 0-terminated dot bracket string
 *
 * Generates all neighbors of structure struc as dot-bracket strings
 * and pushs them on the global stack "stapel".
 */
void RNA2_move_it(char * struc);


/**
 * @brief free data structures of move generator
 */
void RNA2_free(void);


/**
 * @brief Returns pointer to sequence of RNA structure.
 * @return pointer to sequence of RNA structure
 */
/* char const* RNA2_get_seq(); */

/**
 * @brief Get dot-bracket representation of current structure.
 * @return Dot-bracket representation of current structure
 */
/* char const * RNA2_get_form(); */

/**
 * @brief Return length of current sequence.
 * @return Length of current sequence, i.e. the number of bases
 */
/* int RNA2_get_len(); */

#endif // _rnamoves_h
