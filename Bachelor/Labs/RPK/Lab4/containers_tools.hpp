#pragma once

#include <stdexcept>

template<std::size_t C>
struct dims_info
{
	std::size_t data[C];

	static constexpr auto count = C;

	inline std::size_t &operator[](const std::size_t i) { return data[i]; }
	inline const std::size_t &operator[](const std::size_t i) const { return data[i]; }

	inline std::size_t &external_dim() { return data[0]; }
	inline const std::size_t &external_dim() const { return data[0]; }

	inline std::size_t &internal_dim() { return data[C - 1]; }
	inline const std::size_t &internal_dim() const { return data[C - 1]; }
};

template<std::size_t C1, std::size_t C2>
constexpr inline bool operator==(const dims_info<C1> &d1, const dims_info<C2> &d2)
{
	return false;
}

template<std::size_t C>
constexpr inline bool operator==(const dims_info<C> &d1, const dims_info<C> &d2)
{
	for (std::size_t i = 0; i < C; ++i)
		if (d1[i] != d2[i])
			return false;
	return true;
}

template<std::size_t C1, std::size_t C2>
constexpr inline bool operator!=(const dims_info<C1> &d1, const dims_info<C2> &d2)
{
	return true;
}

template<std::size_t C>
constexpr inline bool operator!=(const dims_info<C> &d1, const dims_info<C> &d2)
{
	for (std::size_t i = 0; i < C; ++i)
		if (d1[i] != d2[i])
			return true;
	return false;
}

template<typename container_t>
inline dims_info<1> dims(const container_t &arg)
{
	return { arg.size() };
}

template<typename container_t>
inline container_t init(const dims_info<1> d)
{
	return container_t(d[0]);
}

template<typename result_t, typename arg_t, typename function_t>
inline auto map(const arg_t &args, const function_t f)
{
	auto result = init<result_t>(dims(args));

	auto arg = args.begin();
	auto res = result.begin();
	auto end = result.end();

	while (res != end) *res++ = f(*arg++);

	return result;
}

template<typename result_t, typename fsts_t, typename snds_t, typename function_t>
inline auto map_pair(const fsts_t &fsts, const snds_t &snds, const function_t f)
{
	if (dims(fsts) != dims(snds))
		throw std::invalid_argument("Mismatch of dimensions.");

	auto result = init<result_t>(dims(fsts));

	auto fst = fsts.begin();
	auto snd = snds.begin();
	auto res = result.begin();
	auto end = result.end();

	while (res != end) *res++ = f(*fst++, *snd++);

	return result;
}

template<typename container_t, typename unary_predicate_t>
inline auto find_first_iterator(const container_t &container, const unary_predicate_t p)
{
	auto it = container.begin();
	auto end = container.end();

	while (it != end && !p(*it)) ++it;
	return it;
}