#pragma once

#include <cassert>
#include <vector>
#include <list>

namespace dvs
{
	template <class S>
	class Space;
	template <class S>
	class Vector
	{
		S* _values;
		Space<S>& _space;
	public:
		Vector(S* values, Space<S>& space) noexcept : _values(values), _space(space) { }
		~Vector() noexcept
		{
			if (_values)
				_space._dealloc(_values);
		}
		Vector() = delete;
		Vector(Vector&& b) noexcept : _values(b._values), _space(b._space)
		{
			b._values = nullptr;
		}
		Vector(const Vector& b) noexcept : _values(b._space._alloc()), _space(b._space())
		{
			for (int i = 0; i < _space._dim; ++i)
				_values[i] = b[i];
		}

		int spaceCheck(const Vector& b) const noexcept
		{
			assert(&_space == &b._space);
			return _space._dim;
		}

		Vector& operator=(Vector&& b) noexcept
		{
			assert(&_space == &b._space);
			if (&b == this)
				return *this;
			if (_values)
				_space._dealloc(_values);
			_values = b._values;
			b._values = nullptr;
			return *this;
		}
		Vector& operator=(const Vector& b) noexcept
		{
			if (!_values)
				_values = _space._alloc();
			for (int i = 0, c = spaceCheck(b); i < c; ++i)
				_values[i] = b[i];
			return *this;
		}

		Space<S>& space() const { return _space; }
		int size() const { return _space._dim; }

		S& operator[](int i) { return _values[i]; }
		const S& operator[](int i) const { return _values[i]; }

		S* begin() { return _values; }
		const S* begin() const { return _values; }
		S* end() { return _values + _space._dim; }
		const S* end() const { return _values + _space._dim; }

		Vector& operator+=(const Vector& b)
		{
			for (int i = 0, c = spaceCheck(c); i < c; ++i)
				_values[i] += b[i];
			return *this;
		}
		Vector& operator+=(S b)
		{
			for (int i = 0, c = size(); i < c; ++i)
				_values[i] += b;
			return *this;
		}
		Vector& operator-=(const Vector& b)
		{
			for (int i = 0, c = spaceCheck(b); i < c; ++i)
				_values[i] -= b[i];
			return *this;
		}
		Vector& operator-=(S b)
		{
			for (int i = 0, c = size(); i < c; ++i)
				_values[i] -= b;
			return *this;
		}
		Vector& operator*=(const Vector& b)
		{
			for (int i = 0, c = spaceCheck(b); i < c; ++i)
				_values[i] *= b[i];
			return *this;
		}
		Vector& operator*=(S b)
		{
			for (int i = 0, c = size(); i < c; ++i)
				_values[i] *= b;
			return *this;
		}
		Vector& operator/=(const Vector& b)
		{
			for (int i = 0, c = spaceCheck(b); i < c; ++i)
				_values[i] /= b[i];
			return *this;
		}
		Vector& operator/=(S b)
		{
			for (int i = 0, c = size(); i < c; ++i)
				_values[i] /= b;
			return *this;
		}

		Vector& negate(const Vector& b)
		{
			for (int i = 0, c = spaceCheck(b); i < c; ++i)
				_values[i] = b[i] - _values[i];
			return *this;
		}
		Vector& negate(S b)
		{
			for (int i = 0, c = size(); i < c; ++i)
				_values[i] = b - _values[i];
			return *this;
		}
		Vector& negate()
		{
			for (int i = 0, c = size(); i < c; ++i)
				_values[i] = -_values[i];
			return *this;
		}

		Vector& invert(const Vector& b)
		{
			for (int i = 0, c = spaceCheck(b); i < c; ++i)
				_values[i] = b[i] / _values[i];
			return *this;
		}
		Vector& invert(S b)
		{
			for (int i = 0, c = size(); i < c; ++i)
				_values[i] = b / _values[i];
			return *this;
		}
	};
	template <class S>
	class Space
	{
		static_assert(std::is_trivially_default_constructible_v<S>);
		static_assert(std::is_trivially_destructible_v<S>);
		friend class Vector<S>;

		using Vector = Vector<S>;
	public:
		struct Free
		{
			Free* next;
		};
		struct Block
		{
			Block* prev = nullptr;
			S data[(4096 - sizeof(decltype(prev))) / sizeof(S)];

			Block() = default;
			Block(Block* prev) : prev(prev) { }
		};
		static_assert(sizeof(Block) <= 4096);

		Block _root_block;
		const int _dim;
		const int _blocksize; // The size of a block in number of scalars
		Block* _last_block;
		S* _next = nullptr;
		S* _end = nullptr;
		size_t _in_use = 0;
		Free* _free = nullptr;

		Vector* _zero_basis;

		S* _free_as_S() const { return reinterpret_cast<S*>(_free); }

		S* _advance_next()
		{
			if (_next == _end)
			{
				_last_block = new Block(_last_block);
				_next = _last_block->data;
				_end = _next + _blocksize;
			}
			const auto result = _next;
			_next += _dim;
			return result;
		}

		void _unwind_next()
		{
			while (_free && _free_as_S() + _dim == _next)
			{
				_next = _free_as_S();
				_free = _free->next;
				if (_next + _blocksize == _end && _last_block != &_root_block)
				{
					delete std::exchange(_last_block, _last_block->prev);
					_end = _last_block->data + _blocksize;
					_next = _end;
				}
			}
		}

		S* _alloc()
		{
			++_in_use;
			if (_free)
			{
				const auto result = _free_as_S();
				_free = _free->next;
				_unwind_next();
				return result;
			}
			return _advance_next();
		}
		void _dealloc(S* ptr)
		{
			--_in_use;
			if (ptr + _dim == _next)
			{
				_next = ptr;
				_unwind_next();
			}
			else
			{
				const auto new_free = reinterpret_cast<Free*>(ptr);
				new_free->next = _free;
				_free = new_free;
			}
		}
	public:
		Space(int dim) : _dim(dim), _blocksize((sizeof(Block) / (sizeof(S)*_dim))*_dim), 
			_last_block(&_root_block), _next(_root_block.data), _end(_next + _blocksize)
		{
			assert(sizeof(Free) <= sizeof(S)*_dim);
			_zero_basis = reinterpret_cast<Vector*>(new char[sizeof(Vector)*(_dim + 1)]);
			new (_zero_basis) Vector(_alloc(), *this);
			for (int i = 0; i < _dim; ++i)
				_zero_basis[0][i] = S(0);
			for (int k = 0; k < _dim; ++k)
			{
				new (_zero_basis + k + 1) Vector(_alloc(), *this);
				for (int i = 0; i < _dim; ++i)
					_zero_basis[k+1][i] = i == k ? S(1) : S(0);
			}
		}
		Space() = delete;
		Space(Space&&) = delete;
		Space(const Space&) = delete;
		~Space()
		{
			for (int k = _dim; k > 0; --k)
				_zero_basis[k].~Vector();
			_zero_basis[0].~Vector();
			delete[] reinterpret_cast<char*>(_zero_basis);
			assert(_in_use == 0);
		}
		Space& operator=(Space&&) = delete;
		Space& operator=(const Space&) = delete;

		const Vector& zero() const { return _zero_basis[0]; }
		const Vector& basis(int i) const { return _zero_basis[i + 1]; }

		template <class... Components, class Enable = std::enable_if_t<sizeof...(Components) >= 2>>
		Vector operator()(Components&&... c)
		{
			assert(sizeof...(Components) == _dim);
			Vector result(_alloc(), *this);
			int i = 0;
			((result[i++] = S(std::forward<Components>(c))), ...);
			return result;
		}
		template <class T, class Enable = ifExpr<T>>
		Vector operator()(const T& e)
		{
			assert(this == &e.space());
			return { e };
		}

		int dim() const { return _dim; }
	};

	template <class S> Vector<S> operator-(      Vector<S>&& a) { return std::move(a.negate()); }
	template <class S> Vector<S> operator-(const Vector<S>&& a) { auto r = a; a.negate(); return r; }

	template <class S> Vector<S> operator+(      Vector<S>&& a,       Vector<S>&& b) { return std::move(a += b); }
	template <class S> Vector<S> operator+(      Vector<S>&& a, const Vector<S>&  b) { return std::move(a += b); }
	template <class S> Vector<S> operator+(const Vector<S>&  a,       Vector<S>&& b) { return std::move(b += a); }
	template <class S> Vector<S> operator+(const Vector<S>&  a, const Vector<S>&  b) { auto r = a; r += b; return r; }
	template <class S> Vector<S> operator+(Vector<S>&& a,      S  b) { return std::move(a += b); }
	template <class S> Vector<S> operator+(S  a,      Vector<S>&& b) { return std::move(b += a); }
	template <class S> Vector<S> operator+(const Vector<S>& a, S  b) { auto r = a; r += b; return r; }
	template <class S> Vector<S> operator+(S  a, const Vector<S>& b) { auto r = b; r += a; return r; }

	template <class S> Vector<S> operator-(      Vector<S>&& a,       Vector<S>&& b) { return std::move(a -= b); }
	template <class S> Vector<S> operator-(      Vector<S>&& a, const Vector<S>&  b) { return std::move(a -= b); }
	template <class S> Vector<S> operator-(const Vector<S>&  a,       Vector<S>&& b) { return std::move(b.negate(a)); }
	template <class S> Vector<S> operator-(const Vector<S>&  a, const Vector<S>&  b) { auto r = a; r -= b; return r; }
	template <class S> Vector<S> operator-(Vector<S>&& a,      S  b) { return std::move(a -= b); }
	template <class S> Vector<S> operator-(S  a,      Vector<S>&& b) { return std::move(b.negate(a)); }
	template <class S> Vector<S> operator-(const Vector<S>& a, S  b) { auto r = a; r -= b; return r; }
	template <class S> Vector<S> operator-(S a, const Vector<S>&  b) { auto r = b; r.negate(b); return r; }

	template <class S> Vector<S> operator*(      Vector<S>&& a,       Vector<S>&& b) { return std::move(a *= b); }
	template <class S> Vector<S> operator*(      Vector<S>&& a, const Vector<S>&  b) { return std::move(a *= b); }
	template <class S> Vector<S> operator*(const Vector<S>&  a,       Vector<S>&& b) { return std::move(b *= a); }
	template <class S> Vector<S> operator*(const Vector<S>&  a, const Vector<S>&  b) { auto r = a; r *= b; return r; }
	template <class S> Vector<S> operator*(Vector<S>&& a,      S  b) { return std::move(a *= b); }
	template <class S> Vector<S> operator*(S  a,      Vector<S>&& b) { return std::move(b *= a); }
	template <class S> Vector<S> operator*(const Vector<S>& a, S  b) { auto r = a; r *= b; return r; }
	template <class S> Vector<S> operator*(S  a, const Vector<S>& b) { auto r = b; r *= a; return r; }

	template <class S> Vector<S> operator/(      Vector<S>&& a,       Vector<S>&& b) { return std::move(a /= b); }
	template <class S> Vector<S> operator/(      Vector<S>&& a, const Vector<S>&  b) { return std::move(a /= b); }
	template <class S> Vector<S> operator/(const Vector<S>&  a,       Vector<S>&& b) { return std::move(b.invert(a)); }
	template <class S> Vector<S> operator/(const Vector<S>&  a, const Vector<S>&  b) { auto r = a; r /= b; return r; }
	template <class S> Vector<S> operator/(Vector<S>&& a,      S  b) { return std::move(a /= b); }
	template <class S> Vector<S> operator/(S  a,      Vector<S>&& b) { return std::move(b.invert(a)); }
	template <class S> Vector<S> operator/(const Vector<S>& a, S  b) { auto r = a; r /= b; return r; }
	template <class S> Vector<S> operator/(S a, const Vector<S>&  b) { auto r = b; r.invert(b); return r; }

	template <class S>
	auto dot(const Vector<S>& a, const Vector<S>& b) 
	{
		auto sp = a[0] * b[0];
		for (int i = 1, c = a.spaceCheck(b); i < c; ++i)
			sp += a[i] * b[i];
		return sp;
	}

	template <class S> 
	auto square(const Vector<S>& a)
	{
		auto ss = a[0] * a[0];
		for (int i = 1, c = a.size(); i < c; ++i)
			ss += a[i] * a[i];
		return ss;
	}

	template <class S>
	S length(const Vector<S>& a)
	{
		return sqrt(square(a));
	}
	
	template <class S>
	struct Decomposition
	{
		S length;
		Vector<S> direction;

		Decomposition(Vector<S> a) : direction(std::move(a)) { length = dvs::length(direction); direction /= length; }
	};
	template <class S> Decomposition<S> decompose(      Vector<S>&& a) { return { std::move(a) }; }
	template <class S> Decomposition<S> decompose(const Vector<S>&  a) { return {           a  }; }

	template <class S> Vector<S> direction(      Vector<S>&& a) { return decompose(std::move(a)).direction; }
	template <class S> Vector<S> direction(const Vector<S>&  a) { return decompose(          a ).direction; }
}
