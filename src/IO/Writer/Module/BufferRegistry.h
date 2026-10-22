// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_IO_WRITER_MODULE_BUFFERREGISTRY_H_
#define SEISSOL_SRC_IO_WRITER_MODULE_BUFFERREGISTRY_H_

#include "IO/Writer/Writer.h"

#include <cstddef>
#include <unordered_map>
#include <unordered_set>
#include <vector>

namespace seissol::io::writer::module {

struct BufferPointer {
  std::size_t size{0};
  int id{-1};
};

//! A buffer that the n-th data source of a kind uses in every plan.
struct BufferSlot {
  int id{-1};
  const void* pointer{nullptr};
  std::size_t size{0};
};

/**
 * @brief The buffer operations BufferRegistry needs, i.e. the part of ASYNC it uses.
 *
 * Splitting this out keeps the bookkeeping below independent of ASYNC, which is what makes it
 * testable without a running module.
 */
class BufferAllocator {
  public:
  BufferAllocator() = default;
  virtual ~BufferAllocator() = default;
  BufferAllocator(const BufferAllocator&) = delete;
  BufferAllocator(BufferAllocator&&) = delete;
  auto operator=(const BufferAllocator&) -> BufferAllocator& = delete;
  auto operator=(BufferAllocator&&) -> BufferAllocator& = delete;

  //! @brief Registers a buffer. A null @p pointer asks the allocator to provide the memory.
  virtual int addBuffer(const void* pointer, std::size_t size) = 0;
  virtual void resizeBuffer(int id, const void* pointer, std::size_t size) = 0;
  //! @brief The memory of an allocator-provided buffer; nullptr if it is empty.
  virtual void* managedBuffer(int id) = 0;

  /**
   * @brief Whether all ranks have to add and resize the same buffers, in the same order.
   *
   * Where an executor joins the parts of several ranks under one id (ASYNC mode NODE), adding and
   * resizing a buffer is collective.
   */
  [[nodiscard]] virtual bool collective() const { return false; }

  /**
   * @brief Sets every entry of @p flags that is set on any rank; only needed where collective()
   * holds.
   */
  virtual void agreeOnAny(std::vector<int>& flags) {}
};

/**
 * @brief Decides which buffer each data source of a write plan uses.
 *
 * A plan refers to its data by buffer id, and the executor looks the data up under that id, so
 * this mapping is the whole contract between the two sides of the write. It carries state across
 * writes: a pass-through buffer keeps the id of its pointer, and a managed buffer reuses an
 * existing allocation of the same size rather than adding a new one on every output step.
 */
class BufferRegistry {
  public:
  explicit BufferRegistry(BufferAllocator& allocator) : allocator_(allocator) {}

  /**
   * @brief Assigns an id to every distributed data source of @p writer and returns the ids in
   * plan order, ready to be sent.
   *
   * A source occurring several times in the plan is handled once. Within one call, two sources
   * never share an id, even when they would be interchangeable across calls.
   *
   * Where the allocator is collective, the n-th pass-through source and the n-th managed source of
   * a plan always take the n-th buffer of their kind: the plans of all ranks are alike, so the
   * buffers line up across the ranks, while their pointers and sizes do not. A buffer is resized
   * on all ranks whenever it changes on one of them.
   */
  std::vector<int> assign(Writer& writer);

  private:
  std::vector<int> assignCollective(Writer& writer);
  int idFor(DataSource* source, const std::unordered_set<int>& used);
  int idForPassThrough(WriteBuffer* source, const std::unordered_set<int>& used);
  int idForManaged(AdhocBuffer* source, const std::unordered_set<int>& used);

  BufferAllocator& allocator_;
  //! Pass-through buffers, keyed by the pointer they hand out.
  std::unordered_map<const void*, BufferPointer> pointerMap_;
  //! Managed buffers, keyed by their size, so that an allocation can be reused.
  std::unordered_map<std::size_t, std::vector<int>> bufferMap_;
  //! Where the allocator is collective: the pass-through buffers by their place among their kind.
  std::vector<BufferSlot> passThroughSlots_;
  //! Where the allocator is collective: the managed buffers by their place among their kind.
  std::vector<BufferSlot> managedSlots_;
};

} // namespace seissol::io::writer::module

#endif // SEISSOL_SRC_IO_WRITER_MODULE_BUFFERREGISTRY_H_
