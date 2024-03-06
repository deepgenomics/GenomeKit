/*
Copyright (C) 2016-2023 Deep Genomics Inc. All Rights Reserved.
*/
#include "interval.h"

#include "file.h"
#include "format.h"
#include "gk_assert.h"
#include "strutil.h"
#include <cstring>

BEGIN_NAMESPACE_GK

string interval_t::as_str() const
{
	// TODO: data_dir injected as a context
	return format("{}:{}-{}:{}:{}", get_chrom_names(refg).chrom_as_sv(chrom), start(), end() - 1,
				  strand_as_char(strand), get_refg_registry().refg_as_sv(refg));
}

pos_t as_pos(std::string_view s) { return as_int(s); }

END_NAMESPACE_GK
