#ifndef __COMPAT__
#define __COMPAT__

/*
 * This file provides compatibility functions to make up for differences between compilers.
 *
 * (c) Greg Reynolds 2007, 2008. This file is not GPL and is BSD style licensed.
 */

// If using Microsoft Visual C++.
#ifdef _MSC_VER

// Turn off over-zealous warnings about functions like "sscanf".
#define _CRT_SECURE_NO_WARNING

// These parts of the C++ standard are missing from Microsoft's compilers.
#define or ||
#define and &&
#define not !

#endif // _MSC_VER

namespace compat 
{

//! Finds the minimum of its two arguments.
// This function can't be called "min", even in its own namespace or
// Visual C++ barfs.
template <typename arg_t> arg_t minimum(const arg_t& x, const arg_t& y)
{
	if( x < y )
		return x;
	else
		return y;
}

//! Finds the maximum of its two arguments.
// This function can't be called "max", even in its own namespace or
// Visual C++ barfs.
template <typename arg_t> arg_t maximum(const arg_t& x, const arg_t& y)
{
	if( x > y )
		return x;
	else
		return y;
}

} // namespace compat
#endif
