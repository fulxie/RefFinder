<?php

//
// Edit History:
//
//  $Author: dickmunroe $
//  $Date: 2007/06/16 17:34:33 $
//
//  Dick Munroe (munroe@csworks.com) 04-Jan-2006
//      Initial version created.
//
//  Dick Munroe (munroe@csworks.com) 04-Jan-2006
//      Add random variable generation using the Box-Mueller transformations.
//
//  Dick Munroe (munroe@csworks.com) 05-Jan-2006
//      Add a mean, median, and mode functions.
//
//  Dick Munroe (munroe@csworks.com) 06-Jan-2006
//      Add standard deviation and variance functions.
//
//  Dick Munroe (munroe@csworks.com) 19-Jan-2006
//      Improve the numerical accuracy and performance of the variance and
//      standard deviation functions (suggested by Martin Weis).
//
//  Dick Munroe (munroe@csworks.com) 22-Jan-2006
//      Swap the integration limits internally to the method if they are
//      out of order.
//
//  Dick Munroe (munroe@csworks.com) 25-Jan-2006
//      Add an equation solver using steepest descent and another using
//      bisection.
//
//  Dick Munroe (munroe@csworks.com) 26-Jan-2006
//      Add a limiting feature to the bisection equation solver.
//
//  Dick Munroe (munroe@csworks.com) 27-Jan-2006
//      Add conversion to rational numbers.
//
//  Dick Munroe (munroe@csworks.com) 01-Apr-2007
//      Add find local minimum/maximum.
//
//  Dick Munroe (munroe@csworks.com) 25-May-2007
//      I must have lost the factor function back in the server crash last
//      year.
//
//  Dick Munroe (munroe@csworks.com) 28-May-2007
//      factor negative numbers as well.
//
//  Dick Munroe (munroe@csworks.com) 14-June-2007
//      Add gaussian elimination.
//      Add matrix inversion using gaussian elimination.
//
//  Dick Munroe (munroe@csworks.com) 31-Dec-2010
//      Add functions to do additiona math on arrayes and some simple
//      statistical calculations.
//

/**
 * @author Dick Munroe <munroe@csworks.com>
 * @copyright copyright @ 2006-2010 by Dick Munroe, Cottage Software Works, Inc.
 * @license http://www.csworks.com/publications/ModifiedNetBSD.html
 * @version 1.0.11
 * @package Numerical
 * @example ./example.php
 */

/*
 * Provide a number of numerical analysis and statistical functions.
 */

define("Numerical_Critical_Error", 0, true) ;
define("Numerical_Critical_Minimum", 1, true) ;
define("Numerical_Critical_Maximum", 2, true) ;
define("Numerical_Inversion_Size", 3, true) ;

class Numerical
{
    /**
     * Add two vectors or a constant and a vector producing a vector.
     * Indices missing from either parameter are assumed to be 0.
     * 
     * @access public
     * @param mixed either numeric variable, constant or array.
     * @param mixed either numeric variable, constant or array
     * @return mixed either array or numeric.
     */
    
    function add(
    	$theFirst,
    	$theSecond)
    {
    	if (!is_array($theFirst) && !is_array($theSecond))
    	{
    		return $theFirst + $theSecond ;
    	}
    	
    	if (is_array($theFirst) && !is_array($theSecond))
    	{
    		foreach ($theFirst as $k => $v)
    		{
    			$xxx[$k] = $v + $theSecond ;
    		}
    		
    		return $xxx ;
    	}
    	
    	if (is_array($theSecond) && !is_array($theFirst))
    	{
    		foreach ($theSecond as $k => $v)
    		{
    			$xxx[$k] = $v + $theFirst ;
    		}
    		
    		return $xxx ;
    	}

    	$xxx = array() ;
    	
    	foreach ($theFirst as $k => $v)
    	{
    		if (!isset($theSecond[$k]))
    		{
    			$xxx[$k] = $v ;
    		}
    	}

    	foreach ($theSecond as $k => $v)
    	{
    		if (isset($theFirst[$k]))
    		{
    			$xxx[$k] = $theFirst[$k] + $v ;
    		}
    		else 
    		{
    			$xxx[$k] = $v ;
    		}
    	}
    	
    	return $xxx ;
   }
    
   	/**
   	 * Divide producing a result.
   	 * 
   	 * The arguments can be any mix of array and floating value.  The
   	 * incides of the arrays must match for division to occur.
   	 * Only those indices common to both arrays will produce
   	 * result.  Possible combinations of arguments are:
   	 * 
   	 * 	value / value
   	 * 	array / value
   	 * 	value / array
   	 * 	array / array
   	 * 
   	 * Array indices are preserved.
   	 *
   	 * @param mixde $theFirst Either a number orn an array of numbers.  The dividend.
   	 * @param mixed $theSecond Either a number or an array of numbers.  The divisor.
   	 * @return mixed an array if at least one of the parameters is an array, otherwise
   	 * 				 a floating point value.
   	 */
    function divide(
    	$theFirst,
    	$theSecond)
    {
    	if (!is_array($theFirst) && !is_array($theSecond))
    	{
    		return floatval($theFirst) / floatval($theSecond) ;
    	}

    	if (is_array($theFirst) && !is_array($theSecond))
    	{
    		foreach ($theFirst as $k => $v)
    		{
    			$xxx[$k] = floatval($v) / floatval($theSecond) ;
    		}
    		
    		return $xxx ;
    	}
    	
    	if (is_array($theSecond) && !is_array($theFirst))
    	{
    		foreach ($theSecond as $k => $v)
    		{
    			$xxx[$k] = floatval($theFirst) / floatval($v) ;
    		}
    		
    		return $xxx ;
    	}
    	
			$xxx = array() ;
    	$common = array_intersect_key($theFirst, $theSecond) ;
    	
    	if (!$common)
    	{
    		return $xxx ;
    	}
    	
			foreach ($common as $k => $v)
			{
				$xxx[$k] = floatval($v) / floatval($theSecond[$k]) ;
			}

    	return $xxx ;
   }
    
   	/**
   	 * Calculate the covariance of two data sets.  The indicies of the
   	 * data sets determines the common data between the two data sets.
   	 * The common data is then used to calculate the covariance.
   	 *
   	 * @param array $theXData
   	 * @param array $theYData
   	 * @return float
   	 */
   	
   	function covariance(
   		$theXData,
   		$theYData)
   	{
   		$xData = array_intersect_key($theXData, $theYData) ;
   		$yData = array_intersect_key($theYData, $xData) ;

   		/*
   		 * Either the sample size is exactly 1 in one of the data sets
   		 * or there is no overlap between the data sets.  In either case
   		 * there is no covariance between the data.
   		 */
   		
   		if (count($xData) <= 1)
   		{
   			return 0.0 ;
   		}
   		
   		/*
   		 * See http:/en.wikipedia.org/wiki/Covariance for the details on this
   		 * calculation.
   		 */
   		
   		$xy = Numerical::multiply($xData, $yData) ;
   		
   		$xy = array_sum($xy) / count($xData) ;
   		
   		$x = array_sum($xData) / count($xData) ;
   		
   		$y = array_sum($yData) / count($xData) ;
   		
   		return $xy - ( $x * $y ) ;
   	}
   	
    /**
     * Find the critical point within the specified interval.
     *
     * Critical will not find all critical points within the interval.  The interval must be chosen
     * so that only one critical point exists.  The function must also be continuous within the
     * boundaries.
     *
     * @access public
     * @param float, by reference $theCriticalPoint the value of the critical point
     * @param function $theFunction the function to be evaluated.
     * @param float $theLowLimit The lower limit of the interval
     * @param float $theHighLimit The high limit of the integral
     * @param float [optional] $theDeltaLimit If the interval to the next guess is within this limit,
     *              the next guess is returned as the critical point.
     * @param float [optional] $theDerivativeLimit If the slope of the derivative is within this limit
     *              the current guess is returned as the critical point.
     * @return float Status of the operation (See defines, above).
     */

    function critical(
        &$theCriticalPoint,
        $theFunction,
        $theLowerBound,
        $theUpperBound,
        $theDerivativeInterval = 1e-03,
        $theDeltaLimit = 1e-06,
        $theDerivativeLimit = 1e-06)
    {
        $lhs = $theFunction($theLowerBound) ;
        $xxx = $theFunction($theLowerBound + $theDerivativeInterval) ;
        $lhsDerivative =  ($xxx - $lhs) / $theDerivativeInterval ;

        $rhs = $theFunction($theUpperBound) ;
        $xxx = $theFunction($theUpperBound - $theDerivativeInterval) ;
        $rhsDerivative =  ($xxx - $rhs) / -$theDerivativeInterval ;

        /*
         * A sufficient, but not necessary, condition for the existance of
         * a minimum or maximum within the interval is that the derivative
         * changes sign from one bound to the other.
         */

        if (($lhsDerivative * $rhsDerivatave) > 0)
        {
            return Numerical_Critical_Error ;
        }

        if ($lhsDerivative > 0.0)
        {
            $theCriticalPoint =
                Numerical::criticalMaximum(
                    $theFunction,
                    $theLowerBound,
                    $theUpperBound,
                    $theDerivativeInterval,
                    $theDeltaLimit,
                    $theDerivativeLimit) ;

            return Numerical_Critical_Maximum ;
        }
        else
        {
            $theCriticalPoint =
                Numerical::criticalMinimum(
                    $theFunction,
                    $theLowerBound,
                    $theUpperBound,
                    $theDerivativeInterval,
                    $theDeltaLimit,
                    $theDerivativeLimit) ;

            return Numerical_Critical_Minimum ;
        }
    }

    function criticalMinimum(
        $theFunction,
        $theLowerBound,
        $theUpperBound,
        $theDerivativeInterval,
        $theDeltaLimit,
        $theDerivativeLimit)
    {
        do
        {
            $theDelta = ($theUpperBound - $theLowerBound) / 2.0 ;

            $theGuess = $theLowerBound + $theDelta ;

            if ($theDelta < $theDeltaLimit)
            {
                return $theGuess ;
            }

            $xxx = $theFunction($theGuess + $theDerivativeInterval) ;
            $yyy = $theFunction($theGuess) ;

            $theGuessDerivative = ($xxx - $yyy) / $theDerivativeInterval ;

            if (abs($theGuessDerivative) < $theDerivativeLimit)
            {
                return $theGuess ;
            }

            if ($xxx > $yyy)
            {
                $theUpperBound = $theGuess ;
            }
            else
            {
                $theLowerBound = $theGuess ;
            }
        }
        while (TRUE) ;
    }

    function criticalMaximum(
        $theFunction,
        $theLowerBound,
        $theUpperBound,
        $theDerivativeInterval,
        $theDeltaLimit,
        $theDerivativeLimit)
    {
        do
        {
            $theDelta = ($theUpperBound - $theLowerBound) / 2.0 ;

            $theGuess = $theLowerBound + $theDelta ;

            if ($theDelta < $theDeltaLimit)
            {
                return $theGuess ;
            }

            $xxx = $theFunction($theGuess + $theDerivativeInterval) ;
            $yyy = $theFunction($theGuess) ;

            $theGuessDerivative = ($xxx - $yyy) / $theDerivativeInterval ;

            if (abs($theGuessDerivative) < $theDerivativeLimit)
            {
                return $theGuess ;
            }

            if ($xxx > $yyy)
            {
                $theLowerBound = $theGuess ;
            }
            else
            {
                $theUpperBound = $theGuess ;
            }
        }
        while (TRUE) ;
    }

    /**
     * Produce the prime factors of a number using the sieve of Erastophenes.
     *
     * @desc Factor an integer.
     * @access public
     * @param integer $theNumber the number to be factored.
     * @return array the key is the factor, the value of the array at that point is the number
     *               of times that factor occurs in the output.
     */

    function &factor($theNumber)
    {
        $theNumber = abs($theNumber) ;

        /*
         * The upper bound of the prime factors is the square root of the
         * number.
         */

        $upperBound = ceil(sqrt(floatval($theNumber))) ;

        if ($upperBound == 1)
        {
            $theFactors = array(1 => 1) ;
            return $theFactors ;
        }

        /*
         * Generate a list of all possible factors, excluding 1 which is always
         * a factor.
         */

        for ($i = 2; $i <= $upperBound; $i++)
        {
            $theFactors[$i] = 0 ;
        }

        /*
         * Get possible prime factors.
         */

        for ($i = 2; $i <= $upperBound; $i++)
        {
            if (isset($theFactors[$i]))
            {
                for ($j = $i + $i; $j <= $upperBound; $j = $j + $i)
                {
                    unset($theFactors[$j]) ;
                }
            }
        }

        /*
         * Go through the possible factors, counting the ones
         * that exist, eliminating the ones that don't.
         */

        foreach (array_keys($theFactors) as $aFactor)
        {
            if (($theNumber % $aFactor) == 0)
            {
                do
                {
                    $theNumber = $theNumber / $aFactor ;
                    $theFactors[$aFactor]++ ;
                } while (($theNumber % $aFactor) == 0) ;
            }
            else
            {
                unset($theFactors[$aFactor]) ;
            }
        }

        /*
         * One is always a prime factor and it appears once.
         */

        $theFactors[1] = 1 ;

        /*
         * Any residual amount after all the divisions above is also
         * a prime factor and occurs only once.
         */

        $theFactors[$theNumber] = 1 ;

        /*
         * The original input can be reconstructed by taking the sum
         * of the keys raised to the power of their value.
         */

        return $theFactors ;
    }

    /**
     * Generate a Gaussian normal distribution with a specified mean and standard
     * deviation.  The default distribution generated is the standard
     * normal distribution with a mean of 0.0 and a standard deviation of
     * 1.0.
     *
     * @desc Draw a value from a normal distribution.
     * @access public
     * @param float $x The parameter of the function.
     * @param float $theMean The mean of the distribution, by default 0.0.
     * @param float $theStandardDeviation The standard deviation of the distribution, by default, 1.0.
     * @return float the probability of the x occuring in the distribution.
     */

    function gaussian($x, $theMean = 0.0, $theStandardDeviation = 1.0)
    {
        $part1 = (1/($theStandardDeviation * sqrt(2 * M_PI))) ;
        $part2 = exp(- pow(($x - $theMean), 2) / (2 * pow($theStandardDeviation, 2))) ;
        return $part1 * $part2 ;
    }

    /**
     * Solve a system of equations using gaussian elimination.
     *
     * @access public
     * @param integer Number of unknowns in the equation.
     * @param mixed Row major set of coefficients for the matrix to be solved in augmented form.
     * @return reference to mixed of the array containing the solution.
     */

    function &gaussianElimination($n, $a)
    {
        /*
         * Loop across the columns.  All other loops will work to the diagonal.
         */

        $rr = 0 ;
        $theColumns = count($a[0]) ;

        for ($cc = 0; $cc < $n; $cc++)
        {
            /*
             * Find the column with the largest value.
             */

            $maxIndex = $rr ;

            for ($r = $rr + 1; $r < $n; $r++)
            {
                if (abs($a[$r][$cc]) > abs($a[$maxIndex][$cc]))
                {
                    $maxIndex = $r ;
                }
            }

            /*
             * Swap the row containing the max value for the "current" row.
             */

            if ($rr != $maxIndex)
            {
                $aaa = $a[$rr] ;
                $a[$rr] = $a[$maxIndex] ;
                $a[$maxIndex] = $aaa ;
            }

            /*
             * If scaling would result in a divide by 0, then there is no appropariate scaling
             */

            if ($a[$rr][$cc] != 0.0)
            {
                $aaa =& $a[$rr] ;

                /*
                 * Scale the current row so that it's left will be 1.
                 */

                $factor = $aaa[$rr] ;

                for ($r = $rr; $r < $theColumns; $r++)
                {
                    $aaa[$r] = $aaa[$r]  / $factor ;
                }

                unset($aaa) ;

                /*
                 * For all following rows, subtract the appropriate multiple of the current
                 * maximum row to eliminate the left most non-zero column.
                 */
                for ($r = $rr + 1; $r < $n; $r++)
                {
                    $aaa = $a[$rr] ;

                    $factor = $a[$r][$rr] ;

                    for ($c = $rr; $c < $theColumns; $c++)
                    {
                        $a[$r][$c] = $a[$r][$c] - ( $a[$rr][$c] * $factor ) ;
                    }
                }
            }

            $rr++ ;

            if ($rr == $n)
            {
                break ;
            }
        }

        /*
         * Now apply the solution backwards through all the rows to get to the
         * identity matrix and a solution.
         */

        for ($r = $n - 1; $r > 0; $r--)
        {
            for ($rr = $r - 1; $rr >= 0; $rr--)
            {
                $aaa =& $a[$rr] ;

                $factor = $aaa[$r] ;

                $aaa[$r] = $aaa[$r] - ( $a[$r][$r] * $factor ) ;

                for ($c = $n; $c < $theColumns; $c++)
                {
                    $aaa[$c] = $aaa[$c] - ( $a[$r][$c] * $factor ) ;
                }

                unset($aaa) ;
            }
        }

        return $a ;
    }

    /**
     * @desc invert a matrix using gausian elimination.
     *
     * @access public
     * @param mixed the matrix to be inverted.
     * @return reference to array
     */

    function &gaussianInversion($a)
    {
        $theRows = count($a) ;
        $theColumns = count($a[0]) ;

        if ($theRows != $theColumns)
        {
            return Numerical_Inversion_Size ;
        }

        /*
         * Pad the input array with an identity matrix.
         */

        for ($i = 0; $i < $theRows; $i++)
        {
            $a[$i] = array_pad($a[$i], 2 * $theColumns, 0) ;
            $a[$i][$theColumns + $i] = 1 ;
        }

        $a = Numerical::gaussianElimination($theRows, $a) ;

        for ($i = 0; $i < $theRows; $i++)
        {
            $a[$i] = array_slice($a[$i], $theColumns) ;
        }

        return $a ;
    }

    /**
     * Caluclate the integral of a function over a specified range.
     *
     * Uses the midpoint method of calculation for integral.
     * The function passed in must take one parameter and return the
     * value of the function for that parameter.
     *
     * @access public
     * @param function $theFunction the function to be integrated.
     * @param float $theLowLimit The lower limit of the integral
     * @param float $theHighLimit The high limit of the integral
     * @param integer $theNumberOfSteps The number of steps take to evaluate the integral, defaults to 100.
     * @return float The value of the integral for the specified range.
     */

    function integrate($theFunction, $theLowLimit, $theHighLimit, $theNumberOfSteps = 100)
    {
        if ($theLowLimit > $theHighLimit)
        {
            /*
             * The limits are out of order, swap them.
             */

            $xxx = $theLowLimit ;
            $theLowLimit = $theHighLimit ;
            $theHighLimit = $xxx ;
        }

        $theDelta = ($theHighLimit - $theLowLimit) / $theNumberOfSteps ;
        $theMidpointDelta = $theDelta / 2 ;

        $theArea = 0.0 ;

        for ($i = 0; $i < $theNumberOfSteps; $i++)
        {
            $theValue = $theLowLimit + (($i * $theDelta) + $theMidpointDelta) ;
            $theValue = $theFunction($theValue) ;
            $theArea += $theValue * $theDelta ;
        }

        return $theArea ;
    }
		
    /**
     * Take the natural log of an vector.
     * 
     * Indices are preserved, i.e., the natural log is taken in place.
     *
     * @param mixed $theSamples may be either an array or a single value.
     * @return mixed the natural log(s) of the parameter.
     */
    
    function ln(
    	$theSamples)
  	{
  		if (is_array($theSamples))
  		{
  			foreach ($theSamples as $k => $v)
  			{
  				$theSamples[$k] = log($v) ;
  			}
  		}
  		else 
  		{
  			$theSamples = log($theSamples) ;
  		}
  		
  		return $theSamples ;
  	}

  	/**
     * @desc Calculate the mean of an array of samples.
     * @access public
     * @param array $theSamples An array of samples.
     * @return float the mean of the samples.
     */

    function mean($theSamples)
    {
        return array_sum($theSamples) / floatval(count($theSamples)) ;
    }

    /**
     * @desc Determine the median of a sample.
     * @access public
     * @param array $theSamples An array of samples.
     * @return float the mean of the samples.
     */

    function median($theSamples)
    {
        sort($theSamples) ;

        if ((count($theSamples) % 2) == 0)
        {
            $i = count($theSamples) / 2 ;
            return ($theSamples[$i - 1] + $theSamples[$i]) / 2 ;
        }
        else
        {
            return $theSamples[(int) (count($theSamples) / 2)] ;
        }
    }

    /*
     * @desc Determine the mode[s] of a sample.
     * @access public
     * @param array $theSamples An array of samples.
     * @return mixed the mode[s] of the sample.
     */

    function mode($theSamples)
    {
        $theCounts = array() ;

        $theCount = count($theSamples) ;

        for ($i = 0; $i < $theCount; $i++)
        {
            $theCounts[(string)$theSamples[$i]]++ ;
        }

        arsort($theCounts) ;

        $theModes = array() ;

        foreach ($theCounts as $key => $count)
        {
            if (count($theModes) == 0)
            {
                $theModeCount = $count ;
            }

            if ($theModeCount == $count)
            {
                $theModes[] = $key ;
            }
            else
            {
                break ;
            }
        }

        return $theModes ;
    }

    /**
     * Multiple two vectors together.
     * 
     * Only the indices in common between the two vectors participate
     * in the multiplication.  Multiplication of constant and vector
     * is also supported.
     *
     * @param mixed $theFirst multiplicand.
     * @param mixed $theSecond multiplier.
     * @return mixed The product of the two parameters.
     */
    
		function multiply(
			$theFirst,
			$theSecond)
		{
			if (!is_array($theFirst) && !is_array($theSecond))
			{
				return $theFirst * $theSecond ;
			}
			
			if (is_array($theFirst) && !is_array($theSecond))
			{
				foreach ($theFirst as $k => $v)
				{
					$xxx[$k] = $v * $theSecond ;
				}
				
				return $xxx ;
			}

			if (is_array($theSecond) && !is_array($theFirst))
			{
				foreach ($theSecond as $k => $v)
				{
					$xxx[$k] = $v * $theFirst ;
				}
				
				return $xxx ;
			}
			
			$xxx = array() ;
			$common = array_intersect_key($theFirst, $theSecond) ;
			
			if ($common)
			{
				foreach ($common as $k => $v)
				{
					$xxx[$k] = floatval($v) * floatval($theSecond[$k]) ;
				}
			}
						
			return $xxx ;
		}
		
    /*
     * Note that this algorithm can be unstable for uniformly distributed random numbers
     * close to 0.0.
     *
     * @desc Generate gaussian distributed random numbers using the Box-Mueller basic transform.
     * @param float $theMean The mean of the distribution.
     * @param float $theStandardDeviation The standard deviation of the distribution.
     * @return float A random number.
     */

    function randomGaussianBoxMuellerBasic($theMean = 0.0, $theStandardDeviation = 1.0)
    {
        static $useLast = FALSE ;
        static $y2 ;

        if ($useLast)
        {
            $useLast = FALSE ;
            $y1 = $y2 ;
        }
        else
        {
            /*
             * Generate a uniformly distributed random number in the range
             * (0, 1].
             */

            $theRange = 10000000 ;

            $x1 = mt_rand(1, $theRange) / $theRange ;
            $x2 = mt_rand(1, $theRange) / $theRange ;

            /*
             * Convert this to a pair of gaussian distributed random numbers.
             * The 2nd is returned on the next call.
             */

            $x = sqrt( -2 * log($x1) ) ;

            $y1 = $x * sin( 2 * M_PI * $x2 ) ;
            $y2 = $x * cos( 2 * M_PI * $x2 ) ;

            $useLast = TRUE ;
        }

        return $theMean + $y1 * $theStandardDeviation ;
    }

    /*
     * @desc Generate gaussian distributed random numbers using the Box-Mueller Polar transform.
     * @param float $theMean The mean of the distribution.
     * @param float $theStandardDeviation The standard deviation of the distribution.
     * @return float A random number.
     */

    function randomGaussianBoxMuellerPolar($theMean = 0.0, $theStandardDeviation = 1.0)
    {
        static $useLast = FALSE ;
        static $y2 ;

        if ($useLast)
        {
            $useLast = FALSE ;
            $y1 = $y2 ;
        }
        else
        {
            do
            {
                /*
                 * Generate a pair of uniformly distributed random numbers
                 * in the range [-1..1].
                 */

                $theRange = 10000000 ;

                $x1 = mt_rand(- $theRange, $theRange) / $theRange ;
                $x2 = mt_rand(- $theRange, $theRange) / $theRange ;

                $w = $x1 * $x1 + $x2 * $x2 ;
            } while ($w >= 1.0) ;

            $w = sqrt ( (-2 * log( $w ) ) / $w ) ;
            $y1 = $x1 * $w ;
            $y2 = $x2 * $w ;

            $useLast = TRUE ;
        }

        return $theMean + $y1 * $theStandardDeviation ;
    }

    /*
     * @desc Determine the standard deviation of a sample.
     * @access public
     * @param array $theSamples An array of samples.
     * @return float the standard deviation.
     */

    function standardDeviation($theSamples)
    {
        return sqrt(Numerical::variance($theSamples)) ;
    }

    /*
     * Calculate the rational number version of a floating point number.
     *
     * Note that the rational number returned is NOT guaranteed to be in
     * simplest form, i.e., the algorithm doesn't necessarily find 2/3, but
     * will find 4/6 for it's solution (which is the same as 2/3).
     *
     * A very cool algorithm (and, as near as I can tell, unique).  Basically
     * take a guess, calculate a delta, use the delta to calculate a new
     * denominator, take another guess, keep going until things converge.
     *
     * @param float $x The floating point number to be converted.
     * @param float $theEpsilon How close the fraction should be to x.
     * @param float $theLimit How big the numberator or denominator can get.
     * @return array [0] = numerator [1] = denominator.
     */

    function rational($x, $theEpsilon = 1.0e-06)
    {
        /*
         * Strip off the sign, the algorithm only works on positive
         * numbers.  We'll put the sign back on the way out.
         */

        if ($x < 0.0)
        {
            $theSign = -1 ;
            $x = -$x ;
        }
        else
        {
            $theSign = 1 ;
        }

        /*
         * Strip off the integer portion, we'll add that back in as we return.
         */

        $theInteger = floor($x) ;
        $x = $x - $theInteger ;

        /*
         * Catch the case of a real number with no fractional part.
         */

        if ($x == 0.0)
        {
            return array($theSign * $theInteger, 1) ;
        }

        /*
         * Take an inital guess at the fractional part of the rational number.
         */

        $d = round(1/$x, 0) ;
        $n = round($x * $d, 0) ;

        do
        {
            /*
             * If our current ratio is "close enough", we're done.
             */

            $delta = abs($x - ($n / $d)) ;

            if ($delta < $theEpsilon)
            {
                $n = $theInteger * $d + $n ;

                /*
                 * return error if the numerator overflows.
                 */

                if ($n < 0)
                {
                    return sqrt(-1) ;
                }

                return array($theSign * $n, $d) ;
            }

            /*
             * Figure out the denominator of the fraction representing
             * the delta.
             */

            $d = round(1/$delta, 0) ;

            /*
             * Check to see if we had an integer overflow in calculating
             * the denominator.  If we did, then bail with Nan.
             */

            if ($d < 0)
            {
                return sqrt(-1) ;
            }

            /*
             * The new guess for the numerator is however many pieces of
             * the denominator are necessary.
             */

            $n = round($x * $d, 0) ;
        } while (true) ;
    }

    /*
     * Solve for roots using the Bisection algorithm.
     *
     * The value of the function must change sign from the left to the
     * right hand side.  This method is slow, but certain.
     *
     * @access public
     * @param function $theFunction the function to be solved.
     * @param float $theLHS The left hand side of the bracket.
     * @param float $theRHS The right hand side of the bracket.
     * @param float $theEpsilon How close does the solution have to be to zero to be a solution.
     * @param float $theSSEpsilon How small does the solution space have to be before failure is declared.
     * @return float The value of the input that solves the function.
     */

    function solveBisection($theFunction, $theLHS, $theRHS, $theEpsilon = 1.0e-06, $theSSEpsilon = 1.0e-9)
    {
        while (true)
        {
            /*
             * If the solution space has gotten very small and no
             * solution is within an epsilon of 0.0, then we declare
             * failure.
             */

            if (($theRHS - $theLHS) < $theSSEpsilon)
            {
                return sqrt(-1) ;
            }

            $x = ($theLHS + $theRHS) / 2.0 ;

            /*
             * Calculate the value of the function.
             */

            $y = $theFunction($x) ;

            /*
             * Since this always looks for 0, if we're within epsilon of any 0,
             * quit.
             */

            if (abs($y) < $theEpsilon)
            {
                return $x ;
            }

            /*
             * Switch the appropriate side of the bracket to use
             * the guessed $x value for a root.
             */

            if ($y < 0)
            {
                $theLHS = $x ;
            }
            else
            {
                $theRHS = $x ;
            }
        }
    }

    /*
     * Solve for roots using the steepest descent algorithm.
     *
     * Can iterate a lot for "flat" curves.
     *
     * @access public
     * @param function $theFunction the function to be solved.
     * @param float $theGuess The guess at a solution.
     * @param float $theEpsilon Accuracy of the solution.
     * @param float $theDerivativeDelta Delta over which the derivative is calculated.
     * @return float The value of the input that solves the function.
     */

    function solveSteepestDescent($theFunction, $theGuess, $theEpsilon = 1.0e-06, $theDerivativeDelta = .001,
                                  $theIterationLimit = 20)
    {
        while (true)
        {
            /*
             * Calculate the value of the function.
             */

            $y = $theFunction($theGuess) ;

            /*
             * Since this always looks for 0, if we're within epsilon of any 0,
             * quit.
             */

            if (abs($y) < $theEpsilon)
            {
                return $theGuess ;
            }

            /*
             * There are some function which can take a LONG
             * time to converge, so we limit the number of
             * iterations and return NaN if we get into this
             * situation.
             */

            if ($theIteration++ > $theIterationLimit)
            {
                return @sqrt(-1) ;
            }

            /*
             * The derivative is the slope "at" the point in question ($theGuess).
             * We define the derivative delta as that distance which is considered
             * for all practical purposes, 0.
             */

            $theDerivative = ($theFunction($theGuess + $theDerivativeDelta) - $y) / $theDerivativeDelta ;

            /*
             * The new guess is calculate by assuming that the
             * equation being solved is a line and moving to where
             * the solution would be.  The algegra is pretty easy
             * keeping in mind that the derivative is the slope of
             * the curve at the point $y.  In effect this subtracts
             * the solution for the linear approximation of the
             * curve at the point $theGuess from $theGuess, moving
             * in the direction of the solution.
             */

            $theGuess = $theGuess - ($y / $theDerivative) ;
        }
    }

    /*
     * @desc Determine the variance of a sample.
     * @access public
     * @param array $theSamples An array of samples.
     * @return float the variance.
     */

    function variance($theSamples)
    {
    		if (count($theSamples) == 1)
    		{
    			return 0.0 ;
    		}
    		
        $theMean = Numerical::mean($theSamples) ;
        $theSum = 0.0 ;
        foreach ($theSamples as $theValue)
        {
            $theSum += pow(($theValue - $theMean), 2) ;
        }
        return ($theSum/floatval(count($theSamples) - 1)) ;
     }
}

?>
