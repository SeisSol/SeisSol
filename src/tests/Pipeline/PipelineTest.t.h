#include <cxxtest/TestSuite.h>
#include "Solver/Pipeline/GenericPipeline.h"
#include <string>
#include <memory>
#include <algorithm>
#include <iterator>
#include <iostream>

namespace seissol {
  namespace unit_test {
    class PipelineTest;
  }
}

class seissol::unit_test::PipelineTest : public CxxTest::TestSuite {
private:
  using TestPipeline = seissol::GenericPipeline<4, 1024>;
  struct TestStageCallBack : public TestPipeline::PipelineCallBack {
    explicit TestStageCallBack(std::string& buffer) : buffer(buffer) {}
    void finalize() override {
      buffer.push_back('F');
    }
  protected:
    std::string& buffer;
  };

  struct FirstStage : public TestStageCallBack {
    explicit FirstStage(std::string& buffer) : TestStageCallBack(buffer) {}
    void operator()(size_t, size_t batchSize, size_t) override {
      buffer.push_back('A');
      batchSizes.emplace_back(batchSize);
    }

    std::vector<size_t> getBatchSizes() {return batchSizes;}
  private:
    // Note: it is enough to test batchSizes in only one stage
    std::vector<size_t> batchSizes{};
  };

  struct SecondStage : public TestStageCallBack {
    explicit SecondStage(std::string& buffer) : TestStageCallBack(buffer) {}
    void operator()(size_t, size_t, size_t) override {
      buffer.push_back('B');
    }
  };

  struct ThirdStage : public TestStageCallBack {
    explicit ThirdStage(std::string& buffer) : TestStageCallBack(buffer) {}
    void operator()(size_t, size_t, size_t) override {
      buffer.push_back('C');
    }
  };

  struct FourthStage : public TestStageCallBack {
    explicit FourthStage(std::string& buffer) : TestStageCallBack(buffer) {}
    void operator()(size_t, size_t, size_t) override {
      buffer.push_back('D');
    }
  };

  static std::string format(std::string&& stringWithSpaces) {
    static auto notSpace = [](unsigned char x){return std::isspace(x) == 0;};

    std::string stringWoSpaces{};
    std::copy_if(stringWithSpaces.begin(),
                 stringWithSpaces.end(),
                 std::back_inserter(stringWoSpaces),
                 notSpace);
    return stringWoSpaces;
  }

public:
  void testGeneral() {
    TestPipeline pipeline{};
    TS_ASSERT_EQUALS(TestPipeline::NumStages, 4);
    TS_ASSERT_EQUALS(TestPipeline::TailSize, 3);
  }

  void testSuperShortPipeline() {
    std::string testBuffer{};
    std::array<std::shared_ptr<TestPipeline::PipelineCallBack>, 4> callBacks = {
        std::make_shared<FirstStage>(testBuffer),
        std::make_shared<SecondStage>(testBuffer),
        std::make_shared<ThirdStage>(testBuffer),
        std::make_shared<FourthStage>(testBuffer),
    };
    TestPipeline pipeline{};
    for (unsigned i = 0; i < callBacks.size(); ++i) {
      pipeline.registerCallBack(i, callBacks[i].get());
    }

    pipeline.run(TestPipeline::DefaultBatchSize - 3);

    auto testedBatchSizes = static_cast<FirstStage*>(callBacks[0].get())->getBatchSizes();
    std::vector<size_t> expectedBatchSizes{TestPipeline::DefaultBatchSize - 3};
    TS_ASSERT_EQUALS(testedBatchSizes.size(), expectedBatchSizes.size());
    for (size_t i = 0; i < expectedBatchSizes.size(); ++i)
      TS_ASSERT_EQUALS(testedBatchSizes[i], expectedBatchSizes[i]);

    std::string expectedResults{PipelineTest::format("AF BF CF DF")};
    TS_ASSERT_EQUALS(testBuffer.size(), expectedResults.size());
    TS_ASSERT_SAME_DATA(testBuffer.data(), expectedResults.data(), expectedResults.size());
  }

  void testShortPipeline() {
    std::string testBuffer{};
    std::array<std::shared_ptr<TestPipeline::PipelineCallBack>, 4> callBacks = {
        std::make_shared<FirstStage>(testBuffer),
        std::make_shared<SecondStage>(testBuffer),
        std::make_shared<ThirdStage>(testBuffer),
        std::make_shared<FourthStage>(testBuffer),
    };
    TestPipeline pipeline{};
    for (unsigned i = 0; i < callBacks.size(); ++i) {
      pipeline.registerCallBack(i, callBacks[i].get());
    }

    pipeline.run(2 * TestPipeline::DefaultBatchSize - 3);

    auto testedBatchSizes = static_cast<FirstStage*>(callBacks[0].get())->getBatchSizes();
    std::vector<size_t> expectedBatchSizes{TestPipeline::DefaultBatchSize,
                                           TestPipeline::DefaultBatchSize - 3};
    TS_ASSERT_EQUALS(testedBatchSizes.size(), expectedBatchSizes.size());
    for (size_t i = 0; i < expectedBatchSizes.size(); ++i)
      TS_ASSERT_EQUALS(testedBatchSizes[i], expectedBatchSizes[i]);

    std::string expectedResults{PipelineTest::format("A AFB BFC CFD DF")};
    TS_ASSERT_EQUALS(testBuffer.size(), expectedResults.size());
    TS_ASSERT_SAME_DATA(testBuffer.data(), expectedResults.data(), expectedResults.size());
  }

  void testWithOneFullPipelineIteration() {
    std::string testBuffer{};
    std::array<std::shared_ptr<TestPipeline::PipelineCallBack>, 4> callBacks = {
        std::make_shared<FirstStage>(testBuffer),
        std::make_shared<SecondStage>(testBuffer),
        std::make_shared<ThirdStage>(testBuffer),
        std::make_shared<FourthStage>(testBuffer),
    };
    TestPipeline pipeline{};
    for (unsigned i = 0; i < callBacks.size(); ++i) {
      pipeline.registerCallBack(i, callBacks[i].get());
    }

    pipeline.run(4 * TestPipeline::DefaultBatchSize - 3);

    auto testedBatchSizes = static_cast<FirstStage*>(callBacks[0].get())->getBatchSizes();
    std::vector<size_t> expectedBatchSizes{TestPipeline::DefaultBatchSize,
                                           TestPipeline::DefaultBatchSize,
                                           TestPipeline::DefaultBatchSize,
                                           TestPipeline::DefaultBatchSize - 3};
    TS_ASSERT_EQUALS(testedBatchSizes.size(), expectedBatchSizes.size());
    for (size_t i = 0; i < expectedBatchSizes.size(); ++i)
      TS_ASSERT_EQUALS(testedBatchSizes[i], expectedBatchSizes[i]);

    std::string expectedResults{PipelineTest::format("A AB ABC AFBCD BFCD CFD DF")};
    TS_ASSERT_EQUALS(testBuffer.size(), expectedResults.size());
    TS_ASSERT_SAME_DATA(testBuffer.data(), expectedResults.data(), expectedResults.size());
  }

  void testLongPipeline() {
    std::string testBuffer{};
    std::array<std::shared_ptr<TestPipeline::PipelineCallBack>, 4> callBacks = {
        std::make_shared<FirstStage>(testBuffer),
        std::make_shared<SecondStage>(testBuffer),
        std::make_shared<ThirdStage>(testBuffer),
        std::make_shared<FourthStage>(testBuffer),
    };
    TestPipeline pipeline{};
    for (unsigned i = 0; i < callBacks.size(); ++i) {
      pipeline.registerCallBack(i, callBacks[i].get());
    }

    pipeline.run(4 * TestPipeline::DefaultBatchSize + 3);

    auto testedBatchSizes = static_cast<FirstStage*>(callBacks[0].get())->getBatchSizes();
    std::vector<size_t> expectedBatchSizes{TestPipeline::DefaultBatchSize,
                                           TestPipeline::DefaultBatchSize,
                                           TestPipeline::DefaultBatchSize,
                                           TestPipeline::DefaultBatchSize,
                                           3};
    TS_ASSERT_EQUALS(testedBatchSizes.size(), expectedBatchSizes.size());
    for (size_t i = 0; i < expectedBatchSizes.size(); ++i)
      TS_ASSERT_EQUALS(testedBatchSizes[i], expectedBatchSizes[i]);


    std::string expectedResults{PipelineTest::format("A AB ABC ABCD AFBCD BFCD CFD DF")};
    TS_ASSERT_EQUALS(testBuffer.size(), expectedResults.size());
    TS_ASSERT_SAME_DATA(testBuffer.data(), expectedResults.data(), expectedResults.size());
  }
};
