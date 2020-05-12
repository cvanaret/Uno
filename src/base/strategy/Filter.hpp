#ifndef FILTER_H
#define FILTER_H

#include <ostream>
#include <vector>
#include <list>
#include <map>
#include <memory>

struct FilterConstants {
    double Beta; /*!< Margin around filter */
    double Gamma; /*!< Margin around filter (sloping margin) */
};

struct FilterEntry {
    double infeasibility_measure;
    double optimality_measure;
};

/*! \class Filter
 * \brief Filter
 *
 *  Filter
 */
class Filter {
public:
    Filter(FilterConstants& constants);
    virtual ~Filter();

    std::list<FilterEntry> entries;
    double upper_bound; /*!< Upper bound on constraint violation */
    unsigned int max_size; /*!< Max filter size */
    FilterConstants constants; /*!< Set of constants */

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
    virtual void add(double infeasibility_measure, double optimality_measure);

    /*!
     *  Check whether a point is acceptable
     *  Virtual method (may be overwritten in subclasses)
     * 
     * \param constraint: constraint value
     * \param objective: objective value
     */
    virtual bool accept(double infeasibility_measure, double optimality_measure);

    /*!
     *  Check if a point is acceptable wrt the current point
     *  Virtual method (may be overwritten in subclasses)
     * 
     * \param current_constraint: current constraint value
     * \param current_objective: current objective value
     * \param trial_constraint: trial objective value
     * \param trial_objective: trial objective value
     */
    virtual bool improves_current_iterate(double current_infeasibility_measure, double current_optimality_measure, double trial_infeasibility_measure, double trial_optimality_measure);

    /*!
     *  Compute the actual reduction resulting from taking the step
     *  Virtual method (may be overwritten in subclasses)
     * 
     * \param current_iterate: current iterate and its evaluations
     * \param trial_objective: objective value of the trial point
     */
    virtual double compute_actual_reduction(double current_objective, double current_residual, double trial_objective);

    double eta_min();
    double omega_min();
    
    /*!
     *  Print the filter
     */
    friend std::ostream& operator<<(std::ostream &stream, Filter& filter);
};

/*! \class NonmonotoneFilter
 * \brief Non-monotonic filter
 *
 *  Non-monotonic filter
 */
class NonmonotoneFilter : public Filter {
public:
    NonmonotoneFilter(FilterConstants& constants, int number_dominated_entries = 3);

    /*!
     *  Add a point to the filter
     * 
     * \param constraint: constraint value
     * \param objective: objective value
     */
    void add(double infeasibility_measure, double optimality_measure);

    /*!
     *  Check whether a point is acceptable
     * 
     * \param constraint: constraint value
     * \param objective: objective value
     */
    bool accept(double infeasibility_measure, double optimality_measure);

    /*!
     *  Check if a point is acceptable wrt the current point
     * 
     * \param current_constraint: current constraint value
     * \param current_objective: current objective value
     * \param trial_constraint: trial objective value
     * \param trial_objective: trial objective value
     */
    bool improves_current_iterate(double current_infeasibility_measure, double current_optimality_measure, double trial_infeasibility_measure, double trial_optimality_measure);

    /*!
     *  Compute the actual reduction resulting from taking the step
     * 
     * \param current_iterate: current iterate and its evaluations
     * \param trial_objective: objective value of the trial point
     */
    double compute_actual_reduction(double current_objective, double current_residual, double trial_objective);

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

    int number_dominated_entries; /*!< Memory of filter */
};

class FilterFactory {
public:
    static std::shared_ptr<Filter> create(std::map<std::string, std::string> options);
};

#endif // FILTER_H
