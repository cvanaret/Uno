#ifndef FILTER_H
#define FILTER_H

#include <ostream>
#include <vector>

struct FilterConstants {
	double Beta; /*!< Margin around filter */
	double Gamma; /*!< Margin around filter (sloping margin) */
};

/*! \class Filter
* \brief Filter
*
*  Filter
*/
class Filter {
	public:
		Filter(FilterConstants& constants);
		
		std::vector<double> constraints; /*!< Array of constraint residuals */
		std::vector<double> objective; /*!< Array of objective values */
		double upper_bound; /*!< Upper bound on constraint violation */
		int size; /*!< Current filter size */
		int max_size; /*!< Max filter size */
		FilterConstants constants;/*!< Set of constants */
		
		/*!
         *  Reset filter size to zero
         */
		void reset();

		/*!
         *  Add a point to the filter
         *  Virtual method (may be overwritten in subclasses)
         * 
         * \param constraint: constraint value
         * \param objective: objective value
         */
		virtual void add(double constraint, double objective);
		
		/*!
         *  Check whether a point is acceptable
         *  Virtual method (may be overwritten in subclasses)
         * 
         * \param constraint: constraint value
         * \param objective: objective value
         */
		virtual bool query(double constraint, double objective);

		/*!
         *  Check if a point is acceptable wrt the current point
         *  Virtual method (may be overwritten in subclasses)
         * 
         * \param current_constraint: current constraint value
         * \param current_objective: current objective value
         * \param trial_constraint: trial objective value
         * \param trial_objective: trial objective value
         */
		virtual bool query_current_iterate(double current_constraint, double current_objective, double trial_constraint, double trial_objective);
		
		/*!
         *  Compute the actual reduction resulting from taking the step
         *  Virtual method (may be overwritten in subclasses)
         * 
         * \param current_iterate: current iterate and its evaluations
         * \param trial_objective: objective value of the trial point
         */
		virtual double compute_actual_reduction(double current_objective, double current_residual, double new_objective);

		/*!
         *  Print the filter
         */
		friend std::ostream& operator<< (std::ostream &stream, Filter& filter);
	
	protected:
		/*!
         *  Shift entries left
         * 
         * \param start: start index
         * \param length: length of the slice
         */
		void shift_left(int start, int length);
		
		/*!
         *  Shift entries right
         * 
         * \param start: start index
         * \param length: length of the slice
         */
		void shift_right(int start, int length);
};

/*! \class NonMonotonicFilter
* \brief Non-monotonic filter
*
*  Non-monotonic filter
*/
class NonMonotonicFilter: Filter {
	public:
		NonMonotonicFilter(FilterConstants& constants, int number_dominated_entries = 3);
	
		/*!
         *  Add a point to the filter
         * 
         * \param constraint: constraint value
         * \param objective: objective value
         */  
		void add(double constraint, double objective);    
		
		/*!
         *  Check whether a point is acceptable
         * 
         * \param constraint: constraint value
         * \param objective: objective value
         */
		bool query(double constraint, double objective);
		
		/*!
         *  Check if a point is acceptable wrt the current point
         * 
         * \param current_constraint: current constraint value
         * \param current_objective: current objective value
         * \param trial_constraint: trial objective value
         * \param trial_objective: trial objective value
         */
		bool query_current_iterate(double curc, double curf, double newc, double newf);
		
		/*!
         *  Compute the actual reduction resulting from taking the step
         * 
         * \param current_iterate: current iterate and its evaluations
         * \param trial_objective: objective value of the trial point
         */
		double compute_actual_reduction(double current_objective, double current_residual, double new_objective);
		
		int number_dominated_entries; /*!< Memory of filter */
};

#endif // FILTER_H
