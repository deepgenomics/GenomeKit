/*
Copyright (C) 2016-2023 Deep Genomics Inc. All Rights Reserved.
*/
#include "interval.h"

#include "format.h"
#include <algorithm>
#include <vector>

BEGIN_NAMESPACE_GK

static bool _in_list(interval_t i, const std::vector<interval_t>& list)
{
	// TODO: currently the set of intervals is not indexed, so if we ever want to work with
	// a large allow/exclude, may be worth building an interval_table; doing so would
	// require some changes to the table<T,I> template, since right now table<> assumes it's
	// a view of a memory-mapped file, and therefore does not 'own' its memory.
	return list.end() != std::find_if(begin(list), end(list), [&i](auto l) { return i.overlaps(l); });
}

string interval_t::as_str() const
{
	// TODO: data_dir injected as a context
	return format("{}:{}-{}:{}:{}", get_chrom_names(refg).chrom_as_sv(chrom), start(), end() - 1,
				  strand_as_char(strand), get_refg_registry().refg_as_sv(refg));
}

pos_t as_pos(std::string_view s) { return as_int(s); }

interval_filter::interval_filter(interval_filter::interval_fn_t validate_interval)
: _validate_interval(validate_interval)
{
}

void interval_filter::allow(interval_t i)
{
	if (_validate_interval)
		_validate_interval(i);
	_allow.push_back(i);
}

void interval_filter::exclude(interval_t i)
{
	if (_validate_interval)
		_validate_interval(i);
	_exclude.push_back(i);
}

bool interval_filter::filter(interval_t i) const
{
	return !_in_list(i, _exclude) && (std::empty(_allow) || _in_list(i, _allow));
}

void interval_filter::validate() const
{
	std::for_each(std::begin(_allow), std::end(_allow), _validate_interval);
	std::for_each(std::begin(_exclude), std::end(_exclude), _validate_interval);
}

END_NAMESPACE_GK
