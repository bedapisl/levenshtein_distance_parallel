#ifndef LEVENSHTEIN_IMPLEMENTATION_HPP
#define LEVENSHTEIN_IMPLEMENTATION_HPP

#include <interface.hpp>
#include <exception.hpp>
#include <vector>

#include <omp.h>

#define TASK_SIZE 64

template <typename C, typename S>
struct task_info
{
	task_info() : vertical(TASK_SIZE + 1, 0), horizontal(TASK_SIZE + 1, 0) { }

	std::vector<S> vertical;		//vertical[0] is diagonal
	std::vector<S> horizontal;		//horizontal[0] is diagonal
};

template<typename C = char, typename DIST = std::size_t, bool DEBUG = false>
class EditDistance : public IEditDistance<C, DIST, DEBUG>
{
public:
	
	typedef int_fast32_t S;	
	
	/*
	 * \brief Perform the initialization of the functor (e.g., allocate memory buffers).
	 * \param len1, len2 Lengths of first and second string respectively.
	 */
	virtual void init(DIST len1, DIST len2)
	{
		if(len2 < len1)		//len2 is longer
		{
			std::swap(len1, len2);
			swap_strings = true;
		}	
			
		if((len1 % TASK_SIZE != 0) || (len2 % TASK_SIZE != 0))
			throw std::exception();

		number_of_levels = (len1 / TASK_SIZE) + (len2 / TASK_SIZE) - 1;

		length_in_tasks = len2 / TASK_SIZE;
		height_in_tasks = len1 / TASK_SIZE;
		length_in_chars = len2;
		height_in_chars = len1;

		worst_outcome_cached = len2;
		if(len2 == len1)
		{
			same_length_strings = true;
		}
		else
			same_length_strings = false;

		std::vector<S> initial(TASK_SIZE + 1, 0);
		for(int i=0; i<TASK_SIZE + 1; ++i)
			initial[i] = i;
	
		for(S i=0; i < length_in_tasks + 1; ++i)
		{
			old_carry.push_back(task_info<C, S>());
			new_carry.push_back(task_info<C, S>());
		}

		old_carry[0].vertical = initial;
		old_carry[0].horizontal = initial;

		for(int i=0; i<omp_get_max_threads(); ++i)
		{
			thread_data.push_back(task_info<C, S>());
		}
	}

	/*
	 * \brief Compute the distance between two strings.
	 * \param str1, str2 Strings to be compared.
	 * \result The computed edit distance.
	 */
	virtual DIST compute(const std::vector<C> &str1, const std::vector<C> &str2)	
	{	
		std::vector<C>& left = const_cast<std::vector<C>&>(str2);
		std::vector<C>& down = const_cast<std::vector<C>&>(str1);
		if(swap_strings)
			std::swap(down, left);
		
		S last_task;
		S shift = 0;
		S infinity_columns = 0;
		S infinity_rows = 0;
		bool add_infinity_column = false;
		bool add_infinity_row = false;
		S global_minimum;

		for(S next_level = 1; next_level < number_of_levels + 1; ++next_level)
		{
			global_minimum = infinity;
			if(next_level > height_in_tasks)
			{
				shift = next_level - height_in_tasks;
				if(infinity_columns > 0)
					--infinity_columns;

				if((next_level > length_in_tasks) && (infinity_rows > 0))
					--infinity_rows;
			}

			if(add_infinity_column)
			{
				++infinity_columns;
				add_infinity_column = false;
			}

			if(add_infinity_row)
			{
				++infinity_rows;
				add_infinity_row = false;
			}
			
			last_task = std::min(next_level, length_in_tasks);
			S x, y;
		
			#pragma omp parallel for private(x, y) reduction(min : global_minimum)

			for(S task_number = shift + infinity_columns; task_number < last_task - infinity_rows; ++task_number)
			{
				task_info<C, S> & t = thread_data[omp_get_thread_num()];
				y = (next_level - 1 - task_number) * TASK_SIZE;
				x = (task_number) * TASK_SIZE;
				
				if(!all_too_big(old_carry[task_number], y, x))
				{
					S local_minimum = compute_one_task(old_carry[task_number], t, down, left, y, x);
					if(local_minimum < global_minimum)
						global_minimum = local_minimum;
						
					std::swap(new_carry[task_number].horizontal, t.horizontal);	//this dont have to be done if shift > 0
					std::swap(new_carry[task_number + 1].vertical, t.vertical);
					
					if((task_number == shift + infinity_columns) && (next_level - task_number < height_in_tasks))	
														//prvni task, ale ne posledni radek
					{
						if(task_number == 0)
						{
							for(S i=0; i<TASK_SIZE + 1; ++i)
								new_carry[task_number].vertical[i] = next_level * TASK_SIZE + i;
						}
						else if(infinity_columns > 0)
						{
							for(S i=0; i<TASK_SIZE + 1; ++i)
								new_carry[task_number].vertical[i] = infinity;
						}
					}

					if((task_number == last_task - infinity_rows - 1) && (task_number < length_in_tasks - 1))
														//posledni task, ale ne v poslednim sloupci
					{
						if(infinity_rows > 0)
						{
							for(S i=0; i<TASK_SIZE + 1; ++i)
								new_carry[task_number + 1].horizontal[i] = infinity;
						}
						else
						{
							for(S i=0; i<TASK_SIZE + 1; ++i)
								new_carry[task_number + 1].horizontal[i] = next_level * TASK_SIZE + i;
						}
					}
				}			
				else	
				{
					if(task_number == shift + infinity_columns)	//first task
					{
						add_infinity_column = true;
						for(S i=0; i<TASK_SIZE + 1; ++i)
							new_carry[task_number + 1].vertical[i] = infinity;
					}
					if(task_number == last_task - infinity_rows - 1)		//last task
					{
						add_infinity_row = true;
						for(S i=0; i<TASK_SIZE + 1; ++i)
							new_carry[task_number].horizontal[i] = infinity;
					}
					else
					{
						for(S i=0; i<TASK_SIZE + 1; ++i)
						{
							new_carry[task_number].horizontal[i] = infinity;
							new_carry[task_number + 1].vertical[i] = infinity;
						}
					}
				}
			}
		
			worst_outcome_cached = global_minimum;
			new_carry.swap(old_carry);
		}

		return old_carry[length_in_tasks - 1].horizontal[TASK_SIZE];
	}

private:

	//Vraci kolik bude vysledek, pokud se jiz nic nebude rovnat
	S compute_one_task(task_info<C, S>& in, task_info<C, S>& out, const std::vector<C> &down, const std::vector<C> &left, S y, S x)
	{
		out.vertical[0] = in.horizontal[TASK_SIZE];

		S diagonal;
		S tmp;
		S length;
		S height = 0;
		
		while(height < TASK_SIZE)
		{
			diagonal = in.vertical[height];
			in.horizontal[0] = in.vertical[height + 1];
			length = 0;
			while(length < TASK_SIZE)
			{
				if(down[y + height] == left[x + length])
				{
					++length;
					std::swap(diagonal, in.horizontal[length]);
				}
				else
				{	
					++length;
					tmp = diagonal;
					diagonal = in.horizontal[length];
					in.horizontal[length] = std::min(tmp, std::min(in.horizontal[length - 1], in.horizontal[length])) + 1;
				}
			}
			out.vertical[++height] = in.horizontal[TASK_SIZE];
		}
	
		in.horizontal[0] = in.vertical[TASK_SIZE];
		in.horizontal.swap(out.horizontal);		//writes line to horizontal
		
		S min = infinity;
		
		if(length_in_chars - x > height_in_chars - y)
		{
			for(S i=0; i<TASK_SIZE + 1; ++i)
			{
				min = std::min(min, length_in_chars - x - TASK_SIZE + out.vertical[i]);
			}
			return min;
			
		}
		else if(length_in_chars - x < height_in_chars - y)
		{
			for(S i=0; i<TASK_SIZE + 1; ++i)
			{
				min = std::min(min, height_in_chars - y - TASK_SIZE + out.horizontal[i]);
			}
			return min;
		}
		else
		{
			return out.horizontal[TASK_SIZE] + length_in_chars - x - TASK_SIZE;
		}
	}
	
	bool all_too_big(task_info<C, S>& in, S y, S x)
	{
		if(same_length_strings)
		{
			for(S i=0; i<TASK_SIZE + 1; ++i)
			{
				if((std::abs(x + i - y) + in.horizontal[i] <= worst_outcome())
				|| (std::abs(x - y - i) + in.vertical[i] <= worst_outcome()))
				{
					return false;
				}
			}
		}	
		else
		{
			for(S i=0; i<TASK_SIZE + 1; ++i)
			{
				if((min_heuristic(x + i, y, in.horizontal[i]) <= worst_outcome())
				|| (min_heuristic(x, y + i, in.vertical[i]) <= worst_outcome()))
				{
					return false;
				}
			}
		}
		return true;
	}
	
	/*Vrati kolik nejmene by mohla stat cesta pouzivajici dane policko*/
	S min_heuristic(S x, S y, S value)
	{
		return std::abs(length_in_chars - x - (height_in_chars - y)) + value;
	}

	S worst_outcome()
	{
		return worst_outcome_cached;
	}

	bool swap_strings;
	bool same_length_strings;
	S number_of_levels;
	const S infinity = 999999999;
	S worst_outcome_cached;
	
	S length_in_chars, height_in_chars;
	S height_in_tasks, length_in_tasks;

	std::vector<task_info<C, S>> thread_data;
	std::vector<task_info<C, S>> old_carry;
	std::vector<task_info<C, S>> new_carry;
};


#endif
