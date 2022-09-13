/**
 * @file test_ModelSets.cpp 
 * @author Robin Dowell
 * @brief Testing the distribution classes
 * @date 2022-06-16
 * 
 */
#include "gmock/gmock.h"
#include "ModelSets.h"

/**
 * @brief Note that ModelWrapper will free the model (Bidir or Full) it 
 * is point to when destroyed, therefore it must be allocated through "new"
 */
TEST(ModelWrapper, getMuBidir)
{
  // Arrange: bring SUT to desired state
  Bidirectional *thisModel;
  thisModel = new Bidirectional();
  ModelWrapper sut(thisModel);

  // Assert: Verify the outcome
 EXPECT_EQ(sut.getMu(), thisModel->getMu());
}

TEST(ModelWrapper, getMuFull)
{
  // Arrange: bring SUT to desired state
  FullModel *thisModel;
  thisModel = new FullModel();
  ModelWrapper sut(thisModel);

  // Assert: Verify the outcome
 EXPECT_EQ(sut.getMu(), thisModel->bidir.getMu()); 
}

/** Note that the rest of ModelWrapper are pass throughs,  do these
 * need to be tested?   At moment probably not.  But should have error
 * system for if model is EMPTY, and THAT should be tested.
 * 
 * CURRENTLY THIS IS NOT A TEST!!
 */

TEST(ModelWrapper, modelIsEmpty)
{
  // Arrange: bring SUT to desired state
  perStrandInfo coverage(3.,5.);
  perStrandInfo normalize(10.,30.);
  ModelWrapper sut;

  // Act
  sut.initalizeBounds(2.,1., 0.25, 0.3, 1., 30.);
  sut.resetSufficiencyStats();
  sut.getResponsibility(); 
  sut.getMu();
  sut.updateExpectations(5, coverage, normalize);
  sut.calculateRi('.', coverage);
  sut.updateParameters(30.,2.);

  // Assert: Verify the outcome

  // Should all throw an error?
}


TEST(ModelContainer, initializeWithPriors)
{
  // Arrange
  RawData data;
  data.addDataPoints(10, 14, 4);   
  data.addDataPoints(14, 15, 1);
  data.addDataPoints(8, 10, -2);
  dInterval interval(&data,2,1);
    
  ModelContainer sut(3,BIDIR);

  // Act
  sut.initializeWithPriors(&interval);

  EXPECT_EQ(sut.K, 3.);
}

TEST(ModelContainer, SortByMu)
{
  // Arrange
  ModelContainer sut(3,BIDIR);
  // Act
  sut.SortByMu();
  // Assert
  //Need to test in order.
}

TEST(ModelContainer, resetAllSufficiencyStats)
{
  // Arrange
  ModelContainer sut(3,BIDIR);
  // Act
  sut.resetAllSufficiencyStats();
  // Assert
  // Need to check all reset.
}

TEST(ModelContainer, getAllResponsibilities)
{
  // Arrange
  ModelContainer sut(3,BIDIR);
  // Act
  double ans = sut.getAllResponsibilities();
  // Assert
  EXPECT_EQ(ans, 0.0);
}

TEST(ModelContainer, calculateAllRi)
{
  // Arrange
  perStrandInfo coverage(3.,6.);
  perStrandInfo Ri;
  ModelContainer sut(3,BIDIR);

  // Act
  Ri = sut.calculateAllRi(30,coverage);

  // Assert
  EXPECT_EQ(Ri.forward, 4.);
  EXPECT_EQ(Ri.reverse, 4.);

}

TEST(ModelContainer, updateExpectations)
{
  // Arrange
  perStrandInfo coverage(3.,6.);
  perStrandInfo Ri(300., 24.);
  ModelContainer sut(3,BIDIR);

  //Act
  sut.updateExpectations(30, coverage, Ri);

  // Assert
  // Test that expectations got properly update.
}