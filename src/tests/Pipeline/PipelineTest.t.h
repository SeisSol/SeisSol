// SPDX-FileCopyrightText: 2021-2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause

#include "Solver/Pipeline/GenericPipeline.h"
#include <algorithm>
#include <iostream>
#include <iterator>
#include <memory>
#include <string>

namespace seissol::unit_test {

class PipelineTest {
  public:
  using TestPipeline = seissol::GenericPipeline<4, 1024>;

  struct TestStageCallBack : public TestPipeline::PipelineCallBack {
    explicit TestStageCallBack(std::string& buffer) : buffer(buffer) {}

    void finalize() override { buffer.push_back('F'); }

protected:
    std::string& buffer;
  };

  struct FirstStage : public TestStageCallBack {
    explicit FirstStage(std::string& buffer) : TestStageCallBack(buffer) {}

    void operator()(size_t, size_t batchSize, size_t) override {
      buffer.push_back('A');
      batchSizes.emplace_back(batchSize);
    }

    std::vector<size_t> getBatchSizes() { return batchSizes; }

private:
    // Note: it is enough to test batchSizes in only one stage
    std::vector<size_t> batchSizes{};
  };

  struct SecondStage : public TestStageCallBack {
    explicit SecondStage(std::string& buffer) : TestStageCallBack(buffer) {}

    void operator()(size_t, size_t, size_t) override { buffer.push_back('B'); }
  };

  struct ThirdStage : public TestStageCallBack {
    explicit ThirdStage(std::string& buffer) : TestStageCallBack(buffer) {}

    void operator()(size_t, size_t, size_t) override { buffer.push_back('C'); }
  };

  struct FourthStage : public TestStageCallBack {
    explicit FourthStage(std::string& buffer) : TestStageCallBack(buffer) {}

    void operator()(size_t, size_t, size_t) override { buffer.push_back('D'); }
  };

  static std::string format(std::string&& stringWithSpaces) {
    static auto notSpace = [](unsigned char x) { return std::isspace(x) == 0; };

    std::string stringWoSpaces{};
    std::copy_if(stringWithSpaces.begin(),
                 stringWithSpaces.end(),
                 std::back_inserter(stringWoSpaces),
                 notSpace);
    return stringWoSpaces;
  }
};
TEST_CASE("Sizes of pipeline are correct") {
  REQUIRE(PipelineTest::TestPipeline::NumStages == 4);
  REQUIRE(PipelineTest::TestPipeline::TailSize == 3);
}

TEST_CASE("Ultra short pipeline") {
  std::string testBuffer{};
  std::array<std::shared_ptr<PipelineTest::TestPipeline::PipelineCallBack>, 4> callBacks = {
      std::make_shared<PipelineTest::FirstStage>(testBuffer),
      std::make_shared<PipelineTest::SecondStage>(testBuffer),
      std::make_shared<PipelineTest::ThirdStage>(testBuffer),
      std::make_shared<PipelineTest::FourthStage>(testBuffer),
  };
  PipelineTest::TestPipeline pipeline{};
  for (unsigned i = 0; i < callBacks.size(); ++i) {
    pipeline.registerCallBack(i, callBacks[i].get());
  }

  pipeline.run(PipelineTest::TestPipeline::DefaultBatchSize - 3);

  auto testedBatchSizes =
      static_cast<PipelineTest::FirstStage*>(callBacks[0].get())->getBatchSizes();
  std::vector<size_t> expectedBatchSizes{PipelineTest::TestPipeline::DefaultBatchSize - 3};
  REQUIRE(testedBatchSizes.size() == expectedBatchSizes.size());
  for (size_t i = 0; i < expectedBatchSizes.size(); ++i)
    REQUIRE(testedBatchSizes[i] == expectedBatchSizes[i]);

  std::string expectedResults{PipelineTest::format("AF BF CF DF")};
  REQUIRE(testBuffer.size() == expectedResults.size());
  REQUIRE(testBuffer == expectedResults);
}

TEST_CASE("Short pipeline") {
  std::string testBuffer{};
  std::array<std::shared_ptr<PipelineTest::TestPipeline::PipelineCallBack>, 4> callBacks = {
      std::make_shared<PipelineTest::FirstStage>(testBuffer),
      std::make_shared<PipelineTest::SecondStage>(testBuffer),
      std::make_shared<PipelineTest::ThirdStage>(testBuffer),
      std::make_shared<PipelineTest::FourthStage>(testBuffer),
  };
  PipelineTest::TestPipeline pipeline{};
  for (unsigned i = 0; i < callBacks.size(); ++i) {
    pipeline.registerCallBack(i, callBacks[i].get());
  }

  pipeline.run(2 * PipelineTest::TestPipeline::DefaultBatchSize - 3);

  auto testedBatchSizes =
      static_cast<PipelineTest::FirstStage*>(callBacks[0].get())->getBatchSizes();
  std::vector<size_t> expectedBatchSizes{PipelineTest::TestPipeline::DefaultBatchSize,
                                         PipelineTest::TestPipeline::DefaultBatchSize - 3};
  REQUIRE(testedBatchSizes.size() == expectedBatchSizes.size());
  for (size_t i = 0; i < expectedBatchSizes.size(); ++i)
    REQUIRE(testedBatchSizes[i] == expectedBatchSizes[i]);

  std::string expectedResults{PipelineTest::format("A AFB BFC CFD DF")};
  REQUIRE(testBuffer.size() == expectedResults.size());
  REQUIRE(testBuffer == expectedResults);
}

TEST_CASE("One full pipeline iteration") {
  std::string testBuffer{};
  std::array<std::shared_ptr<PipelineTest::TestPipeline::PipelineCallBack>, 4> callBacks = {
      std::make_shared<PipelineTest::FirstStage>(testBuffer),
      std::make_shared<PipelineTest::SecondStage>(testBuffer),
      std::make_shared<PipelineTest::ThirdStage>(testBuffer),
      std::make_shared<PipelineTest::FourthStage>(testBuffer),
  };
  PipelineTest::TestPipeline pipeline{};
  for (unsigned i = 0; i < callBacks.size(); ++i) {
    pipeline.registerCallBack(i, callBacks[i].get());
  }

  pipeline.run(4 * PipelineTest::TestPipeline::DefaultBatchSize - 3);

  auto testedBatchSizes =
      static_cast<PipelineTest::FirstStage*>(callBacks[0].get())->getBatchSizes();
  std::vector<size_t> expectedBatchSizes{PipelineTest::TestPipeline::DefaultBatchSize,
                                         PipelineTest::TestPipeline::DefaultBatchSize,
                                         PipelineTest::TestPipeline::DefaultBatchSize,
                                         PipelineTest::TestPipeline::DefaultBatchSize - 3};
  REQUIRE(testedBatchSizes.size() == expectedBatchSizes.size());
  REQUIRE(testedBatchSizes == expectedBatchSizes);

  std::string expectedResults{PipelineTest::format("A AB ABC AFBCD BFCD CFD DF")};
  REQUIRE(testBuffer.size() == expectedResults.size());
  REQUIRE(testBuffer == expectedResults);
}

TEST_CASE("Long pipeline") {
  std::string testBuffer{};
  std::array<std::shared_ptr<PipelineTest::TestPipeline::PipelineCallBack>, 4> callBacks = {
      std::make_shared<PipelineTest::FirstStage>(testBuffer),
      std::make_shared<PipelineTest::SecondStage>(testBuffer),
      std::make_shared<PipelineTest::ThirdStage>(testBuffer),
      std::make_shared<PipelineTest::FourthStage>(testBuffer),
  };
  PipelineTest::TestPipeline pipeline{};
  for (unsigned i = 0; i < callBacks.size(); ++i) {
    pipeline.registerCallBack(i, callBacks[i].get());
  }

  pipeline.run(4 * PipelineTest::TestPipeline::DefaultBatchSize + 3);

  auto testedBatchSizes =
      static_cast<PipelineTest::FirstStage*>(callBacks[0].get())->getBatchSizes();
  std::vector<size_t> expectedBatchSizes{PipelineTest::TestPipeline::DefaultBatchSize,
                                         PipelineTest::TestPipeline::DefaultBatchSize,
                                         PipelineTest::TestPipeline::DefaultBatchSize,
                                         PipelineTest::TestPipeline::DefaultBatchSize,
                                         3};
  REQUIRE(testedBatchSizes.size() == expectedBatchSizes.size());
  REQUIRE(testedBatchSizes == expectedBatchSizes);

  std::string expectedResults{PipelineTest::format("A AB ABC ABCD AFBCD BFCD CFD DF")};
  REQUIRE(testBuffer.size() == expectedResults.size());
  REQUIRE(testBuffer == expectedResults);
}

} // namespace seissol::unit_test
