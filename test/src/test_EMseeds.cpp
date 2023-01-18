/**
 * @file test_EMseeds.cpp 
 * @author Robin Dowell
 * @brief Testing the distribution classes
 * @date 2022-06-16
 * 
 */
#include "gmock/gmock.h"
#include "split.h"

#include "EMseeds.h"
#include "Data.h"

TEST(Seeds, writeEmpty) {
  // Arrange
  Seeds sut;
  // Act
  std::string output = sut.writeSeedsAsBedFields();
  //Assert
  EXPECT_EQ(output, "0.\t\t");

}

TEST(Seeds, EmptySize) {
  // Arrange
  Seeds sut;
  // Act
  //Assert
  EXPECT_EQ(sut.getNumSeeds(), 0);
}

class SeedsTest: public :: testing::Test {
    protected:
        void SetUp() override {
            PointCov seed1((double)3., (double).13);
            PointCov seed2((double)18., (double).23);
            PointCov seed3((double)25., (double).03);
            sut.mu_seeds.push_back(seed1);
            sut.mu_seeds.push_back(seed2);
            sut.mu_seeds.push_back(seed3);
        }
    // void TearDown() override { }
    Seeds sut;
};

TEST_F(SeedsTest, getNumSeeds)
{
    // Act: call methods on SUT, capture output
    int numS = sut.getNumSeeds();

    // Assert: Verify the outcome
    EXPECT_EQ(numS, 3);
}

TEST_F(SeedsTest, getMaxSeed)
{
    // Act: call methods on SUT, capture output
    double max = sut.getMaxWeight();

    // Assert: Verify the outcome
    EXPECT_EQ(max, 0.23);
}

TEST_F(SeedsTest, getMinSeed)
{
    // Act: call methods on SUT, capture output
    double min = sut.getMinWeight();

    // Assert: Verify the outcome
    EXPECT_EQ(min, .03);
}

TEST_F(SeedsTest, sortByWeight)
{
    // Act: call methods on SUT, capture output
    sut.SortByWeights();

    // Assert: Verify the outcome
    EXPECT_EQ(sut.mu_seeds[0].coverage, .03);;
}

TEST_F(SeedsTest, sortByCoordinate)
{
    // Act: call methods on SUT, capture output
    sut.SortByPositions();

    // Assert: Verify the outcome
    EXPECT_EQ(sut.mu_seeds[0].coordinate, 3.);;
}

TEST_F(SeedsTest, writeAsBedFields)
{
    // Act: call methods on SUT, capture output
    std::string output = sut.writeSeedsAsBedFields();

    // Assert: Verify the outcome
    EXPECT_EQ(output, "3.\t13.,23.,3.\t3.,18.,25.");

}

TEST_F(SeedsTest, IOonSeeds)
{
    std::string output = sut.writeSeedsAsBedFields();

    std::vector<std::string> lineArray; // Contents of line, split on tab (\t)
    lineArray = string_split(output, '\t');
    std::string seedNum = lineArray[0];
    std::string seedWeights = lineArray[1];
    std::string seedPos = lineArray[2];

    Seeds newSUT;

    // Act: call methods on SUT, capture output
    newSUT.getSeedsfromBedFields(seedNum, seedWeights, seedPos);

    // Assert: Verify the outcome
    EXPECT_EQ(newSUT.getNumSeeds(), 3);
    EXPECT_EQ(newSUT.getMaxWeight(), 0.23);
}

TEST(SeedManager, emptyStringInit)
{
   // Arrange
   SeedManager manager;
   // Act
   std::string output = manager.write_out();
   // Assert
   EXPECT_EQ(output, "");
}

TEST(SeedManager, NulldInterval)
{
   // Arrange
   SeedManager manager;
   manager.setSeeds = new Seeds();
   // Act
   std::string output = manager.write_out();
   // Assert
   EXPECT_EQ(output, "\nSeeds: None");
}

TEST(SeedManager, setupRandomSeeds) 
{
  // Arrange
  SeedManager manager;

  // Act
  std::vector<PointCov> ranseeds = manager.setupRandomSeeds(5, 100);
  std::string seedoutput = tfit::write_VectorPointCov(ranseeds);

  // Assert
  EXPECT_EQ(ranseeds.size(), 5);
  // Note these are "fixed" output because of test seeding is 1974
  EXPECT_EQ(seedoutput, "[46.1,1.0] [76.5,1.0] [22.8,1.0] [74.6,1.0] [17.6,1.0] ");
}

TEST(SeedManager, weightRandomly) 
{
  // Arrange
  SeedManager manager;
  std::vector<PointCov> ranseeds = manager.setupRandomSeeds(5, 100);

  // Act
  manager.weightRandomly(&ranseeds);
  std::string seedoutput = tfit::write_VectorPointCov(ranseeds);

  // Assert
  EXPECT_EQ(ranseeds.size(), 5);
  EXPECT_EQ(seedoutput, "[46.1,0.2] [76.5,0.1] [22.8,0.8] [74.6,0.8] [17.6,0.9] ");
}

TEST(SeedManager, shuffleSeeds) 
{
  // Arrange
  SeedManager manager;
  std::vector<PointCov> ranseeds = manager.setupRandomSeeds(10, 100);
  manager.weightRandomly(&ranseeds);

  manager.setSeeds = new Seeds;
  manager.setSeeds->mu_seeds = ranseeds;

  // Act
  manager.shuffleSeeds();

  std::string postShuffle = tfit::write_VectorPointCov(manager.setSeeds->mu_seeds);
  std::string expectedOut = "[17.6,0.7] [74.6,0.5] [89.0,0.3] [81.1,0.6] [46.1,0.1] [25.8,0.0] [22.8,0.5] [76.5,0.9] [12.1,0.4] [99.5,0.4] ";

  // Assert
  EXPECT_EQ(ranseeds.size(), 10);
  EXPECT_EQ(postShuffle, expectedOut);
}

TEST(SeedManager, addUncertainty) 
{
  // Arrange
  SeedManager manager;
  std::vector<PointCov> ranseeds = manager.setupRandomSeeds(1, 10);
  manager.setSeeds = new Seeds;
  manager.setSeeds->mu_seeds = ranseeds;

  // Act
  double noisymu = manager.addUncertainty(&ranseeds[0],3);

  std::string beforeMu = tfit::prettyDecimal(ranseeds[0].coordinate,1);
  std::string afterMu = tfit::prettyDecimal(noisymu,1);

  // Assert
  EXPECT_EQ(beforeMu, "4.6");
  EXPECT_EQ(afterMu, "2.3");
}


class SeedManagerNoInterval: public :: testing::Test {
  protected:
  void SetUp() override {
    data.addDataPoints(4, 14, 4);   
    data.addDataPoints(14, 23, 3);
    data.addDataPoints(4, 14, -2);   
  }
  RawData data;
};


TEST_F(SeedManagerNoInterval, setupDataLink_NogInterval) 
{
  // Arrange
  dInterval transformedData(&data,2,1);
  SeedManager manager;

  // Act
  manager.setupDataLink(&transformedData);

  std::string output = manager.write_out();
  std::string correctout = 
   "\nData: NO bed4!:min:4.00:max:23.00\tFS: 19.00\tRS: 10.00\tParams bins:10., delta:2.0000, scale:1.000\nSeeds: None";

  // Assert
  EXPECT_EQ(output, correctout);

}

class SeedManagerWithInterval: public :: testing::Test {
  protected:
  void SetUp() override {
    data.addDataPoints(4, 14, 4);   
    data.addDataPoints(14, 23, 3);
    data.addDataPoints(4, 14, -2);   
    data.addDataPoints(14, 23, -6);
    /*
    temp.data = &data;
    data.belongsTo = &temp;
    */
  }
  RawData data;
  bed12 temp = bed12("TestName", 100, 1000, "chrTest", 30, ".", 
                "100\t1000\t0,0,0\t5\t22,30,45,10,5\t50,500,600,800,900"); 
};

TEST_F(SeedManagerWithInterval, setupDataLink_withBed12) 
{
  // Arrange
  dInterval transformedData(&data,2,1);
  SeedManager manager;

  // Act
  manager.setupDataLink(&transformedData);

  std::string output = manager.write_out();
  std::string correctout = "\nData: chrTest:min:4.00:max:23.00\tFS: 19.00\tRS: 19.00\tParams bins:10., delta:2.0000, scale:1.000\nSeeds: 5.\t22.,30.,45.,10.,5.\t50.,500.,600.,800.,900.";

  // Assert
  EXPECT_EQ(output, correctout);
}

TEST_F(SeedManagerWithInterval, setupDataLink_withNoSeeds) 
{
  // Arrange
  temp.seeds = NULL;        // remove previous seeds 
  dInterval transformedData(&data,2,1);
  SeedManager manager;

  // Act
  manager.setupDataLink(&transformedData);

  std::string output = manager.write_out();
  std::string correctout = "\nData: chrTest:min:4.00:max:23.00\tFS: 19.00\tRS: 19.00\tParams bins:10., delta:2.0000, scale:1.000\nSeeds: None";

  // Assert
  EXPECT_EQ(output, correctout);
}

// std::vector<PointCov> grabSeedSet(int K);

TEST_F(SeedManagerWithInterval, grabSeedSetKgreater) 
{
  // Arrange
  dInterval transformedData(&data,2,1);
  SeedManager manager;
  manager.setupDataLink(&transformedData);
  transformedData.bins = 1000;      // Testing within range.

  // Act
  std::vector<PointCov> moreSeeds = manager.grabSeedSet(8);

  std::string newSeedoutput = tfit::write_VectorPointCov(moreSeeds);
  std::string expectedout = "[900.0,0.0] [461.0,1.0] [600.0,0.4] [800.0,0.1] [50.0,0.2] [500.0,0.3] [228.8,1.0] [765.2,1.0] ";

  // Assert
  EXPECT_EQ(moreSeeds.size(), 8); 
  EXPECT_EQ(newSeedoutput, expectedout); 
}


