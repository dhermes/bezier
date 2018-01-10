! Licensed under the Apache License, Version 2.0 (the "License");
! you may not use this file except in compliance with the License.
! You may obtain a copy of the License at
!
!     https://www.apache.org/licenses/LICENSE-2.0
!
! Unless required by applicable law or agreed to in writing, software
! distributed under the License is distributed on an "AS IS" BASIS,
! WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
! See the License for the specific language governing permissions and
! limitations under the License.

module functional_curve

  use, intrinsic :: iso_c_binding, only: c_double
  use types, only: dp
  use functional_test_helpers, only: intersect_and_check
  implicit none
  private &
       case1, &
       case2, &
       case3, &
       case4, &
       case5, &
       case6, &
       case7, &
       case8, &
       case9, &
       case10, &
       case11, &
       case12, &
       case13, &
       case14, &
       case15, &
       case16, &
       case17, &
       case18, &
       case19, &
       case20, &
       case21, &
       case22, &
       case23, &
       case24, &
       case25, &
       case26, &
       case27, &
       case28, &
       case29, &
       case30, &
       case31, &
       case32, &
       case33, &
       case34, &
       case35, &
       case36, &
       case37, &
       case38, &
       case39, &
       case40, &
       case41, &
       case42
  public all_cases

contains

  subroutine case1()

    real(c_double) :: nodes1(3, 2)
    real(c_double) :: nodes2(3, 2)
    real(c_double) :: expected_params(2, 2)

    nodes1(1, :) = [0.0_dp, 0.0_dp]
    nodes1(2, :) = [0.5_dp, 1.0_dp]
    nodes1(3, :) = [1.0_dp, 0.0_dp]

    nodes2(1, :) = [1.125_dp, 0.5_dp]
    nodes2(2, :) = [0.625_dp, -0.5_dp]
    nodes2(3, :) = [0.125_dp, 0.5_dp]

    expected_params(1, :) = [0.21451472732312363_dp, 0.9104852726768764_dp]
    expected_params(2, :) = [0.9104852726768764_dp, 0.21451472732312363_dp]

    call intersect_and_check( &
         1, &
         3, nodes1, &
         3, nodes2, &
         2, expected_params)

  end subroutine case1

  subroutine case2()

    real(c_double) :: nodes3(3, 2)
    real(c_double) :: nodes4(3, 2)
    real(c_double) :: expected_params(2, 2)

    nodes3(1, :) = [0.0_dp, 0.0_dp]
    nodes3(2, :) = [1.5_dp, 3.0_dp]
    nodes3(3, :) = [3.0_dp, 0.0_dp]

    nodes4(1, :) = [3.0_dp, 1.5_dp]
    nodes4(2, :) = [2.625_dp, -0.90625_dp]
    nodes4(3, :) = [-0.75_dp, 2.4375_dp]

    expected_params(1, :) = [0.25_dp, 0.875_dp]
    expected_params(2, :) = [0.75_dp, 0.25_dp]

    call intersect_and_check( &
         2, &
         3, nodes3, &
         3, nodes4, &
         2, expected_params)

  end subroutine case2

  subroutine case3()

    real(c_double) :: nodes1(3, 2)
    real(c_double) :: nodes5(3, 2)
    real(c_double) :: expected_params(2, 2)

    nodes1(1, :) = [0.0_dp, 0.0_dp]
    nodes1(2, :) = [0.5_dp, 1.0_dp]
    nodes1(3, :) = [1.0_dp, 0.0_dp]

    nodes5(1, :) = [0.0_dp, 0.75_dp]
    nodes5(2, :) = [0.5_dp, -0.25_dp]
    nodes5(3, :) = [1.0_dp, 0.75_dp]

    expected_params(1, :) = [0.25_dp, 0.75_dp]
    expected_params(2, :) = [0.25_dp, 0.75_dp]

    call intersect_and_check( &
         3, &
         3, nodes1, &
         3, nodes5, &
         2, expected_params)

  end subroutine case3

  subroutine case4()

    real(c_double) :: nodes1(3, 2)
    real(c_double) :: nodes6(3, 2)
    real(c_double) :: expected_params(2, 1)

    nodes1(1, :) = [0.0_dp, 0.0_dp]
    nodes1(2, :) = [0.5_dp, 1.0_dp]
    nodes1(3, :) = [1.0_dp, 0.0_dp]

    nodes6(1, :) = [0.0_dp, 1.0_dp]
    nodes6(2, :) = [0.5_dp, 0.0_dp]
    nodes6(3, :) = [1.0_dp, 1.0_dp]

    expected_params(1, :) = [0.5_dp]
    expected_params(2, :) = [0.5_dp]

    call intersect_and_check( &
         4, &
         3, nodes1, &
         3, nodes6, &
         1, expected_params)

  end subroutine case4

  subroutine case5()

    real(c_double) :: nodes1(3, 2)
    real(c_double) :: nodes7(3, 2)
    real(c_double) :: expected_params(2, 2)

    nodes1(1, :) = [0.0_dp, 0.0_dp]
    nodes1(2, :) = [0.5_dp, 1.0_dp]
    nodes1(3, :) = [1.0_dp, 0.0_dp]

    nodes7(1, :) = [0.0_dp, 0.265625_dp]
    nodes7(2, :) = [0.5_dp, 0.234375_dp]
    nodes7(3, :) = [1.0_dp, 0.265625_dp]

    expected_params(1, :) = [0.15184468808860432_dp, 0.8481553119113957_dp]
    expected_params(2, :) = [0.15184468808860432_dp, 0.8481553119113957_dp]

    call intersect_and_check( &
         5, &
         3, nodes1, &
         3, nodes7, &
         2, expected_params)

  end subroutine case5

  subroutine case6()

    real(c_double) :: nodes1(3, 2)
    real(c_double) :: nodes8(2, 2)
    real(c_double) :: expected_params(2, 2)

    nodes1(1, :) = [0.0_dp, 0.0_dp]
    nodes1(2, :) = [0.5_dp, 1.0_dp]
    nodes1(3, :) = [1.0_dp, 0.0_dp]

    nodes8(1, :) = [0.0_dp, 0.375_dp]
    nodes8(2, :) = [1.0_dp, 0.375_dp]

    expected_params(1, :) = [0.25_dp, 0.75_dp]
    expected_params(2, :) = [0.25_dp, 0.75_dp]

    call intersect_and_check( &
         6, &
         3, nodes1, &
         2, nodes8, &
         2, expected_params)

  end subroutine case6

  subroutine case7()

    real(c_double) :: nodes1(3, 2)
    real(c_double) :: nodes9(2, 2)
    real(c_double) :: expected_params(2, 1)

    nodes1(1, :) = [0.0_dp, 0.0_dp]
    nodes1(2, :) = [0.5_dp, 1.0_dp]
    nodes1(3, :) = [1.0_dp, 0.0_dp]

    nodes9(1, :) = [0.5_dp, 0.0_dp]
    nodes9(2, :) = [0.5_dp, 0.75_dp]

    expected_params(1, :) = [0.5_dp]
    expected_params(2, :) = [0.6666666666666666_dp]

    call intersect_and_check( &
         7, &
         3, nodes1, &
         2, nodes9, &
         1, expected_params)

  end subroutine case7

  subroutine case8()

    real(c_double) :: nodes10(3, 2)
    real(c_double) :: nodes11(2, 2)
    real(c_double) :: expected_params(2, 1)

    nodes10(1, :) = [0.0_dp, 0.0_dp]
    nodes10(2, :) = [4.5_dp, 9.0_dp]
    nodes10(3, :) = [9.0_dp, 0.0_dp]

    nodes11(1, :) = [0.0_dp, 8.0_dp]
    nodes11(2, :) = [6.0_dp, 0.0_dp]

    expected_params(1, :) = [0.3333333333333333_dp]
    expected_params(2, :) = [0.5_dp]

    call intersect_and_check( &
         8, &
         3, nodes10, &
         2, nodes11, &
         1, expected_params)

  end subroutine case8

  subroutine case9()

    real(c_double) :: nodes8(2, 2)
    real(c_double) :: nodes9(2, 2)
    real(c_double) :: expected_params(2, 1)

    nodes8(1, :) = [0.0_dp, 0.375_dp]
    nodes8(2, :) = [1.0_dp, 0.375_dp]

    nodes9(1, :) = [0.5_dp, 0.0_dp]
    nodes9(2, :) = [0.5_dp, 0.75_dp]

    expected_params(1, :) = [0.5_dp]
    expected_params(2, :) = [0.5_dp]

    call intersect_and_check( &
         9, &
         2, nodes8, &
         2, nodes9, &
         1, expected_params)

  end subroutine case9

  subroutine case10()

    real(c_double) :: nodes1(3, 2)
    real(c_double) :: nodes13(5, 2)
    real(c_double) :: expected_params(2, 4)

    nodes1(1, :) = [0.0_dp, 0.0_dp]
    nodes1(2, :) = [0.5_dp, 1.0_dp]
    nodes1(3, :) = [1.0_dp, 0.0_dp]

    nodes13(1, :) = [0.0_dp, 0.0_dp]
    nodes13(2, :) = [0.25_dp, 2.0_dp]
    nodes13(3, :) = [0.5_dp, -2.0_dp]
    nodes13(4, :) = [0.75_dp, 2.0_dp]
    nodes13(5, :) = [1.0_dp, 0.0_dp]

    expected_params(1, :) = [0.0_dp, 0.3110177634953864_dp, 0.6889822365046137_dp, 1.0_dp]
    expected_params(2, :) = [0.0_dp, 0.3110177634953864_dp, 0.6889822365046137_dp, 1.0_dp]

    call intersect_and_check( &
         10, &
         3, nodes1, &
         5, nodes13, &
         4, expected_params)

  end subroutine case10

  subroutine case11()

    real(c_double) :: nodes14(3, 2)
    real(c_double) :: nodes15(3, 2)
    real(c_double) :: expected_params(2, 1)

    nodes14(1, :) = [0.0_dp, 0.0_dp]
    nodes14(2, :) = [0.375_dp, 0.75_dp]
    nodes14(3, :) = [0.75_dp, 0.375_dp]

    nodes15(1, :) = [0.25_dp, 0.625_dp]
    nodes15(2, :) = [0.625_dp, 0.25_dp]
    nodes15(3, :) = [1.0_dp, 1.0_dp]

    expected_params(1, :) = [0.6666666666666666_dp]
    expected_params(2, :) = [0.3333333333333333_dp]

    call intersect_and_check( &
         11, &
         3, nodes14, &
         3, nodes15, &
         1, expected_params)

  end subroutine case11

  subroutine case12()

    real(c_double) :: nodes14(3, 2)
    real(c_double) :: nodes16(3, 2)
    real(c_double) :: expected_params(2, 2)

    nodes14(1, :) = [0.0_dp, 0.0_dp]
    nodes14(2, :) = [0.375_dp, 0.75_dp]
    nodes14(3, :) = [0.75_dp, 0.375_dp]

    nodes16(1, :) = [0.25_dp, 0.5625_dp]
    nodes16(2, :) = [0.625_dp, 0.1875_dp]
    nodes16(3, :) = [1.0_dp, 0.9375_dp]

    expected_params(1, :) = [0.5_dp, 0.8333333333333334_dp]
    expected_params(2, :) = [0.16666666666666666_dp, 0.5_dp]

    call intersect_and_check( &
         12, &
         3, nodes14, &
         3, nodes16, &
         2, expected_params)

  end subroutine case12

  subroutine case13()

    real(c_double) :: nodes10(3, 2)
    real(c_double) :: nodes17(3, 2)
    real(c_double) :: expected_params(2, 1)

    nodes10(1, :) = [0.0_dp, 0.0_dp]
    nodes10(2, :) = [4.5_dp, 9.0_dp]
    nodes10(3, :) = [9.0_dp, 0.0_dp]

    nodes17(1, :) = [11.0_dp, 8.0_dp]
    nodes17(2, :) = [7.0_dp, 10.0_dp]
    nodes17(3, :) = [3.0_dp, 4.0_dp]

    expected_params(1, :) = [0.3333333333333333_dp]
    expected_params(2, :) = [1.0_dp]

    call intersect_and_check( &
         13, &
         3, nodes10, &
         3, nodes17, &
         1, expected_params)

  end subroutine case13

  subroutine case14()

    real(c_double) :: nodes1(3, 2)
    real(c_double) :: nodes18(3, 2)
    real(c_double) :: expected_params(2, 1)

    nodes1(1, :) = [0.0_dp, 0.0_dp]
    nodes1(2, :) = [0.5_dp, 1.0_dp]
    nodes1(3, :) = [1.0_dp, 0.0_dp]

    nodes18(1, :) = [1.0_dp, 0.0_dp]
    nodes18(2, :) = [1.5_dp, -1.0_dp]
    nodes18(3, :) = [2.0_dp, 0.0_dp]

    expected_params(1, :) = [1.0_dp]
    expected_params(2, :) = [0.0_dp]

    call intersect_and_check( &
         14, &
         3, nodes1, &
         3, nodes18, &
         1, expected_params)

  end subroutine case14

  subroutine case15()

    real(c_double) :: nodes1(3, 2)
    real(c_double) :: nodes19(3, 2)
    real(c_double) :: expected_params(2, 1)

    nodes1(1, :) = [0.0_dp, 0.0_dp]
    nodes1(2, :) = [0.5_dp, 1.0_dp]
    nodes1(3, :) = [1.0_dp, 0.0_dp]

    nodes19(1, :) = [2.0_dp, 0.0_dp]
    nodes19(2, :) = [1.5_dp, 1.0_dp]
    nodes19(3, :) = [1.0_dp, 0.0_dp]

    expected_params(1, :) = [1.0_dp]
    expected_params(2, :) = [1.0_dp]

    call intersect_and_check( &
         15, &
         3, nodes1, &
         3, nodes19, &
         1, expected_params)

  end subroutine case15

  subroutine case16()

    real(c_double) :: nodes1(3, 2)
    real(c_double) :: nodes20(3, 2)
    real(c_double) :: expected_params(2, 4)

    nodes1(1, :) = [0.0_dp, 0.0_dp]
    nodes1(2, :) = [0.5_dp, 1.0_dp]
    nodes1(3, :) = [1.0_dp, 0.0_dp]

    nodes20(1, :) = [1.0_dp, 0.0_dp]
    nodes20(2, :) = [-1.0_dp, 0.25_dp]
    nodes20(3, :) = [1.0_dp, 0.5_dp]

    expected_params(1, :) = [0.09549150281252629_dp, 0.25_dp, 0.6545084971874737_dp, 1.0_dp]
    expected_params(2, :) = [0.3454915028125263_dp, 0.75_dp, 0.9045084971874737_dp, 0.0_dp]

    call intersect_and_check( &
         16, &
         3, nodes1, &
         3, nodes20, &
         4, expected_params)

  end subroutine case16

  subroutine case17()

    real(c_double) :: nodes20(3, 2)
    real(c_double) :: nodes21(3, 2)
    real(c_double) :: expected_params(2, 4)

    nodes20(1, :) = [1.0_dp, 0.0_dp]
    nodes20(2, :) = [-1.0_dp, 0.25_dp]
    nodes20(3, :) = [1.0_dp, 0.5_dp]

    nodes21(1, :) = [-0.125_dp, -0.28125_dp]
    nodes21(2, :) = [0.5_dp, 1.28125_dp]
    nodes21(3, :) = [1.125_dp, -0.28125_dp]

    expected_params(1, :) = [0.0_dp, 0.3454915028125263_dp, 0.75_dp, 0.9045084971874737_dp]
    expected_params(2, :) = [0.9_dp, 0.17639320225002103_dp, 0.3_dp, 0.623606797749979_dp]

    call intersect_and_check( &
         17, &
         3, nodes20, &
         3, nodes21, &
         4, expected_params)

  end subroutine case17

  subroutine case18()

    real(c_double) :: nodes21(3, 2)
    real(c_double) :: nodes22(3, 2)
    real(c_double) :: expected_params(2, 4)

    nodes21(1, :) = [-0.125_dp, -0.28125_dp]
    nodes21(2, :) = [0.5_dp, 1.28125_dp]
    nodes21(3, :) = [1.125_dp, -0.28125_dp]

    nodes22(1, :) = [1.5625_dp, -0.0625_dp]
    nodes22(2, :) = [-1.5625_dp, 0.25_dp]
    nodes22(3, :) = [1.5625_dp, 0.5625_dp]

    expected_params(1, :) = [0.17639320225002103_dp, 0.3_dp, 0.623606797749979_dp, 0.9_dp]
    expected_params(2, :) = [0.37639320225002104_dp, 0.7_dp, 0.823606797749979_dp, 0.1_dp]

    call intersect_and_check( &
         18, &
         3, nodes21, &
         3, nodes22, &
         4, expected_params)

  end subroutine case18

  subroutine case19()

    real(c_double) :: nodes10(3, 2)
    real(c_double) :: nodes23(2, 2)
    real(c_double) :: expected_params(2, 1)

    nodes10(1, :) = [0.0_dp, 0.0_dp]
    nodes10(2, :) = [4.5_dp, 9.0_dp]
    nodes10(3, :) = [9.0_dp, 0.0_dp]

    nodes23(1, :) = [3.0_dp, 4.5_dp]
    nodes23(2, :) = [8.0_dp, 4.5_dp]

    expected_params(1, :) = [0.5_dp]
    expected_params(2, :) = [0.3_dp]

    call intersect_and_check( &
         19, &
         3, nodes10, &
         2, nodes23, &
         1, expected_params)

  end subroutine case19

  subroutine case20()

    real(c_double) :: nodes1(3, 2)
    real(c_double) :: nodes24(3, 2)
    real(c_double) :: expected_params(2, 2)

    nodes1(1, :) = [0.0_dp, 0.0_dp]
    nodes1(2, :) = [0.5_dp, 1.0_dp]
    nodes1(3, :) = [1.0_dp, 0.0_dp]

    nodes24(1, :) = [0.25_dp, 0.375_dp]
    nodes24(2, :) = [0.75_dp, 0.875_dp]
    nodes24(3, :) = [1.25_dp, -0.625_dp]

    expected_params(1, :) = [0.25_dp, 1.0_dp]
    expected_params(2, :) = [0.0_dp, 0.75_dp]

    call intersect_and_check( &
         20, &
         3, nodes1, &
         3, nodes24, &
         2, expected_params)

  end subroutine case20

  subroutine case21()

    real(c_double) :: nodes15(3, 2)
    real(c_double) :: nodes25(4, 2)
    real(c_double) :: expected_params(2, 1)

    nodes15(1, :) = [0.25_dp, 0.625_dp]
    nodes15(2, :) = [0.625_dp, 0.25_dp]
    nodes15(3, :) = [1.0_dp, 1.0_dp]

    nodes25(1, :) = [0.0_dp, 0.5_dp]
    nodes25(2, :) = [0.25_dp, 1.0_dp]
    nodes25(3, :) = [0.75_dp, 1.5_dp]
    nodes25(4, :) = [1.0_dp, 0.5_dp]

    expected_params(1, :) = [0.8587897065534015_dp]
    expected_params(2, :) = [0.873452854150878_dp]

    call intersect_and_check( &
         21, &
         3, nodes15, &
         4, nodes25, &
         1, expected_params)

  end subroutine case21

  subroutine case22()

    real(c_double) :: nodes11(2, 2)
    real(c_double) :: nodes26(4, 2)
    real(c_double) :: expected_params(2, 3)

    nodes11(1, :) = [0.0_dp, 8.0_dp]
    nodes11(2, :) = [6.0_dp, 0.0_dp]

    nodes26(1, :) = [0.375_dp, 7.0_dp]
    nodes26(2, :) = [2.125_dp, 8.0_dp]
    nodes26(3, :) = [3.875_dp, 0.0_dp]
    nodes26(4, :) = [5.625_dp, 1.0_dp]

    expected_params(1, :) = [0.11416126713641388_dp, 0.5_dp, 0.8858387328635862_dp]
    expected_params(2, :) = [0.059041448155901566_dp, 0.5_dp, 0.9409585518440984_dp]

    call intersect_and_check( &
         22, &
         2, nodes11, &
         4, nodes26, &
         3, expected_params)

  end subroutine case22

  subroutine case23()

    real(c_double) :: nodes8(2, 2)
    real(c_double) :: nodes27(4, 2)
    real(c_double) :: expected_params(2, 2)

    nodes8(1, :) = [0.0_dp, 0.375_dp]
    nodes8(2, :) = [1.0_dp, 0.375_dp]

    nodes27(1, :) = [0.125_dp, 0.25_dp]
    nodes27(2, :) = [0.375_dp, 0.75_dp]
    nodes27(3, :) = [0.625_dp, 0.0_dp]
    nodes27(4, :) = [0.875_dp, 0.1875_dp]

    expected_params(1, :) = [0.20998193826534545_dp, 0.4482997212940395_dp]
    expected_params(2, :) = [0.11330925102046059_dp, 0.4310662950587193_dp]

    call intersect_and_check( &
         23, &
         2, nodes8, &
         4, nodes27, &
         2, expected_params)

  end subroutine case23

  subroutine case24()

    real(c_double) :: nodes28(3, 2)
    real(c_double) :: nodes29(3, 2)
    real(c_double) :: expected_params(2, 1)

    nodes28(1, :) = [0.0_dp, 0.0_dp]
    nodes28(2, :) = [-0.5_dp, 1.5_dp]
    nodes28(3, :) = [1.0_dp, 1.0_dp]

    nodes29(1, :) = [-1.0_dp, 1.0_dp]
    nodes29(2, :) = [0.5_dp, 0.5_dp]
    nodes29(3, :) = [0.0_dp, 2.0_dp]

    expected_params(1, :) = [0.5_dp]
    expected_params(2, :) = [0.5_dp]

    call intersect_and_check( &
         24, &
         3, nodes28, &
         3, nodes29, &
         1, expected_params)

  end subroutine case24

  subroutine case25()

    real(c_double) :: nodes29(3, 2)
    real(c_double) :: nodes30(2, 2)
    real(c_double) :: expected_params(2, 1)

    nodes29(1, :) = [-1.0_dp, 1.0_dp]
    nodes29(2, :) = [0.5_dp, 0.5_dp]
    nodes29(3, :) = [0.0_dp, 2.0_dp]

    nodes30(1, :) = [0.5_dp, 0.5_dp]
    nodes30(2, :) = [-0.25_dp, 1.25_dp]

    expected_params(1, :) = [0.5_dp]
    expected_params(2, :) = [0.6666666666666666_dp]

    call intersect_and_check( &
         25, &
         3, nodes29, &
         2, nodes30, &
         1, expected_params)

  end subroutine case25

  subroutine case26()

    real(c_double) :: nodes8(2, 2)
    real(c_double) :: nodes23(2, 2)
    real(c_double) :: expected_params(2, 0)

    nodes8(1, :) = [0.0_dp, 0.375_dp]
    nodes8(2, :) = [1.0_dp, 0.375_dp]

    nodes23(1, :) = [3.0_dp, 4.5_dp]
    nodes23(2, :) = [8.0_dp, 4.5_dp]

    call intersect_and_check( &
         26, &
         2, nodes8, &
         2, nodes23, &
         0, expected_params)

  end subroutine case26

  subroutine case27()

    real(c_double) :: nodes11(2, 2)
    real(c_double) :: nodes31(2, 2)
    real(c_double) :: expected_params(2, 0)

    nodes11(1, :) = [0.0_dp, 8.0_dp]
    nodes11(2, :) = [6.0_dp, 0.0_dp]

    nodes31(1, :) = [0.25_dp, 16.75_dp]
    nodes31(2, :) = [12.25_dp, 0.75_dp]

    call intersect_and_check( &
         27, &
         2, nodes11, &
         2, nodes31, &
         0, expected_params)

  end subroutine case27

  subroutine case28()

    real(c_double) :: nodes32(3, 2)
    real(c_double) :: nodes33(3, 2)
    real(c_double) :: expected_params(2, 1)

    nodes32(1, :) = [-0.25_dp, -0.25_dp]
    nodes32(2, :) = [0.1875_dp, -0.25_dp]
    nodes32(3, :) = [0.625_dp, -0.25_dp]

    nodes33(1, :) = [-0.125_dp, 0.0_dp]
    nodes33(2, :) = [0.0625_dp, -0.5_dp]
    nodes33(3, :) = [0.0_dp, -1.0_dp]

    expected_params(1, :) = [0.23214285714285715_dp]
    expected_params(2, :) = [0.25_dp]

    call intersect_and_check( &
         28, &
         3, nodes32, &
         3, nodes33, &
         1, expected_params)

  end subroutine case28

  subroutine case29()

    real(c_double) :: nodes34(2, 2)
    real(c_double) :: nodes35(2, 2)
    real(c_double) :: expected_params(2, 1)

    nodes34(1, :) = [0.5_dp, 1.875_dp]
    nodes34(2, :) = [0.5_dp, 2.625_dp]

    nodes35(1, :) = [0.625_dp, 2.625_dp]
    nodes35(2, :) = [-0.625_dp, 2.625_dp]

    expected_params(1, :) = [1.0_dp]
    expected_params(2, :) = [0.1_dp]

    call intersect_and_check( &
         29, &
         2, nodes34, &
         2, nodes35, &
         1, expected_params)

  end subroutine case29

  subroutine case30()

    real(c_double) :: nodes36(2, 2)
    real(c_double) :: nodes37(2, 2)
    real(c_double) :: expected_params(2, 1)

    nodes36(1, :) = [0.5_dp, -0.375_dp]
    nodes36(2, :) = [0.0_dp, -0.75_dp]

    nodes37(1, :) = [-0.5_dp, -0.875_dp]
    nodes37(2, :) = [1.25_dp, 0.0_dp]

    expected_params(1, :) = [0.0_dp]
    expected_params(2, :) = [0.5714285714285714_dp]

    call intersect_and_check( &
         30, &
         2, nodes36, &
         2, nodes37, &
         1, expected_params)

  end subroutine case30

  subroutine case31()

    real(c_double) :: nodes38(3, 2)
    real(c_double) :: nodes39(3, 2)
    real(c_double) :: expected_params(2, 1)

    nodes38(1, :) = [0.5_dp, 0.125_dp]
    nodes38(2, :) = [1.25_dp, -0.25_dp]
    nodes38(3, :) = [2.0_dp, 0.5_dp]

    nodes39(1, :) = [0.5_dp, -0.125_dp]
    nodes39(2, :) = [1.0_dp, 0.125_dp]
    nodes39(3, :) = [1.5_dp, -0.125_dp]

    expected_params(1, :) = [0.3333333333333333_dp]
    expected_params(2, :) = [0.5_dp]

    call intersect_and_check( &
         31, &
         3, nodes38, &
         3, nodes39, &
         1, expected_params)

  end subroutine case31

  subroutine case32()

    real(c_double) :: nodes40(2, 2)
    real(c_double) :: nodes41(2, 2)
    real(c_double) :: expected_params(2, 0)

    nodes40(1, :) = [-0.28040801468046_dp, -0.21145489702592957_dp]
    nodes40(2, :) = [-0.30573520623357997_dp, -0.07183588647159315_dp]

    nodes41(1, :) = [-0.2867761039288504_dp, -0.17635312979393736_dp]
    nodes41(2, :) = [-0.3727388074986012_dp, -0.2689285277745333_dp]

    call intersect_and_check( &
         32, &
         2, nodes40, &
         2, nodes41, &
         0, expected_params)

  end subroutine case32

  subroutine case33()

    real(c_double) :: nodes42(4, 2)
    real(c_double) :: nodes43(4, 2)
    real(c_double) :: expected_params(2, 2)

    nodes42(1, :) = [0.0_dp, 2.0_dp]
    nodes42(2, :) = [-0.5_dp, 1.0_dp]
    nodes42(3, :) = [-0.25_dp, 0.75_dp]
    nodes42(4, :) = [-0.09375_dp, 0.828125_dp]

    nodes43(1, :) = [-0.09375_dp, 0.828125_dp]
    nodes43(2, :) = [0.0625_dp, 0.90625_dp]
    nodes43(3, :) = [0.125_dp, 1.3125_dp]
    nodes43(4, :) = [-0.75_dp, 1.625_dp]

    expected_params(1, :) = [0.2546440075000701_dp, 1.0_dp]
    expected_params(2, :) = [0.7453559924999299_dp, 0.0_dp]

    call intersect_and_check( &
         33, &
         4, nodes42, &
         4, nodes43, &
         2, expected_params)

  end subroutine case33

  subroutine case34()

    real(c_double) :: nodes44(4, 2)
    real(c_double) :: nodes45(4, 2)
    real(c_double) :: expected_params(2, 1)

    nodes44(1, :) = [0.0_dp, 2.0_dp]
    nodes44(2, :) = [-0.25_dp, 1.5_dp]
    nodes44(3, :) = [-0.3125_dp, 1.1875_dp]
    nodes44(4, :) = [-0.29296875_dp, 1.009765625_dp]

    nodes45(1, :) = [-0.29296875_dp, 1.009765625_dp]
    nodes45(2, :) = [-0.2734375_dp, 0.83203125_dp]
    nodes45(3, :) = [-0.171875_dp, 0.7890625_dp]
    nodes45(4, :) = [-0.09375_dp, 0.828125_dp]

    expected_params(1, :) = [1.0_dp]
    expected_params(2, :) = [0.0_dp]

    call intersect_and_check( &
         34, &
         4, nodes44, &
         4, nodes45, &
         1, expected_params)

  end subroutine case34

  subroutine case35()

    real(c_double) :: nodes46(4, 2)
    real(c_double) :: nodes47(4, 2)
    real(c_double) :: expected_params(2, 1)

    nodes46(1, :) = [-0.09375_dp, 0.828125_dp]
    nodes46(2, :) = [-0.015625_dp, 0.8671875_dp]
    nodes46(3, :) = [0.0390625_dp, 0.98828125_dp]
    nodes46(4, :) = [-0.03515625_dp, 1.138671875_dp]

    nodes47(1, :) = [-0.03515625_dp, 1.138671875_dp]
    nodes47(2, :) = [-0.109375_dp, 1.2890625_dp]
    nodes47(3, :) = [-0.3125_dp, 1.46875_dp]
    nodes47(4, :) = [-0.75_dp, 1.625_dp]

    expected_params(1, :) = [1.0_dp]
    expected_params(2, :) = [0.0_dp]

    call intersect_and_check( &
         35, &
         4, nodes46, &
         4, nodes47, &
         1, expected_params)

  end subroutine case35

  subroutine case36()

    real(c_double) :: nodes48(2, 2)
    real(c_double) :: nodes49(2, 2)
    real(c_double) :: expected_params(2, 1)

    nodes48(1, :) = [1.0_dp, -3.0_dp]
    nodes48(2, :) = [1.0_dp, 3.0_dp]

    nodes49(1, :) = [10.0_dp, 0.0_dp]
    nodes49(2, :) = [0.0_dp, 0.0_dp]

    expected_params(1, :) = [0.5_dp]
    expected_params(2, :) = [0.9_dp]

    call intersect_and_check( &
         36, &
         2, nodes48, &
         2, nodes49, &
         1, expected_params)

  end subroutine case36

  subroutine case37()

    real(c_double) :: nodes50(2, 2)
    real(c_double) :: nodes54(3, 2)
    real(c_double) :: expected_params(2, 1)

    nodes50(1, :) = [2.25_dp, 3.03125_dp]
    nodes50(2, :) = [2.125_dp, 2.875_dp]

    nodes54(1, :) = [2.1015625_dp, 2.875_dp]
    nodes54(2, :) = [2.16015625_dp, 2.875_dp]
    nodes54(3, :) = [2.220703125_dp, 2.875_dp]

    expected_params(1, :) = [1.0_dp]
    expected_params(2, :) = [0.1993377410829988_dp]

    call intersect_and_check( &
         37, &
         2, nodes50, &
         3, nodes54, &
         1, expected_params)

  end subroutine case37

  subroutine case38()

    real(c_double) :: nodes51(2, 2)
    real(c_double) :: nodes54(3, 2)
    real(c_double) :: expected_params(2, 1)

    nodes51(1, :) = [2.125_dp, 2.875_dp]
    nodes51(2, :) = [2.3125_dp, 2.84375_dp]

    nodes54(1, :) = [2.1015625_dp, 2.875_dp]
    nodes54(2, :) = [2.16015625_dp, 2.875_dp]
    nodes54(3, :) = [2.220703125_dp, 2.875_dp]

    expected_params(1, :) = [0.0_dp]
    expected_params(2, :) = [0.1993377410829988_dp]

    call intersect_and_check( &
         38, &
         2, nodes51, &
         3, nodes54, &
         1, expected_params)

  end subroutine case38

  subroutine case39()

    real(c_double) :: nodes52(3, 2)
    real(c_double) :: nodes54(3, 2)
    real(c_double) :: expected_params(2, 1)

    nodes52(1, :) = [2.25_dp, 3.03125_dp]
    nodes52(2, :) = [2.1875_dp, 2.953125_dp]
    nodes52(3, :) = [2.125_dp, 2.875_dp]

    nodes54(1, :) = [2.1015625_dp, 2.875_dp]
    nodes54(2, :) = [2.16015625_dp, 2.875_dp]
    nodes54(3, :) = [2.220703125_dp, 2.875_dp]

    expected_params(1, :) = [1.0_dp]
    expected_params(2, :) = [0.1993377410829988_dp]

    call intersect_and_check( &
         39, &
         3, nodes52, &
         3, nodes54, &
         1, expected_params)

  end subroutine case39

  subroutine case40()

    real(c_double) :: nodes53(3, 2)
    real(c_double) :: nodes54(3, 2)
    real(c_double) :: expected_params(2, 1)

    nodes53(1, :) = [2.125_dp, 2.875_dp]
    nodes53(2, :) = [2.21875_dp, 2.859375_dp]
    nodes53(3, :) = [2.3125_dp, 2.84375_dp]

    nodes54(1, :) = [2.1015625_dp, 2.875_dp]
    nodes54(2, :) = [2.16015625_dp, 2.875_dp]
    nodes54(3, :) = [2.220703125_dp, 2.875_dp]

    expected_params(1, :) = [0.0_dp]
    expected_params(2, :) = [0.1993377410829988_dp]

    call intersect_and_check( &
         40, &
         3, nodes53, &
         3, nodes54, &
         1, expected_params)

  end subroutine case40

  subroutine case41()

    real(c_double) :: nodes1(3, 2)
    real(c_double) :: nodes55(2, 2)
    real(c_double) :: expected_params(2, 1)

    nodes1(1, :) = [0.0_dp, 0.0_dp]
    nodes1(2, :) = [0.5_dp, 1.0_dp]
    nodes1(3, :) = [1.0_dp, 0.0_dp]

    nodes55(1, :) = [0.0_dp, 0.5_dp]
    nodes55(2, :) = [1.0_dp, 0.5_dp]

    expected_params(1, :) = [0.5_dp]
    expected_params(2, :) = [0.5_dp]

    call intersect_and_check( &
         41, &
         3, nodes1, &
         2, nodes55, &
         1, expected_params)

  end subroutine case41

  subroutine case42()

    real(c_double) :: nodes56(3, 2)
    real(c_double) :: nodes57(3, 2)
    real(c_double) :: expected_params(2, 1)

    nodes56(1, :) = [12.0_dp, 4.0_dp]
    nodes56(2, :) = [-4.0_dp, -4.0_dp]
    nodes56(3, :) = [-4.0_dp, 4.0_dp]

    nodes57(1, :) = [6.0_dp, 1.0_dp]
    nodes57(2, :) = [-2.0_dp, -1.0_dp]
    nodes57(3, :) = [-2.0_dp, 1.0_dp]

    expected_params(1, :) = [0.5_dp]
    expected_params(2, :) = [0.5_dp]

    call intersect_and_check( &
         42, &
         3, nodes56, &
         3, nodes57, &
         1, expected_params)

  end subroutine case42

  subroutine all_cases()

    call case1()
    call case2()
    call case3()
    call case4()
    call case5()
    call case6()
    call case7()
    call case8()
    call case9()
    call case10()
    call case11()
    call case12()
    call case13()
    call case14()
    call case15()
    call case16()
    call case17()
    call case18()
    call case19()
    call case20()
    call case21()
    call case22()
    call case23()
    call case24()
    call case25()
    call case26()
    call case27()
    call case28()
    call case29()
    call case30()
    call case31()
    call case32()
    call case33()
    call case34()
    call case35()
    call case36()
    call case37()
    call case38()
    call case39()
    call case40()
    call case41()
    call case42()

  end subroutine all_cases

end module functional_curve
