#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
pip経由でインストールしたpairlistをテストするスクリプト
"""

import numpy as np

def test_import():
    """基本的なインポートテスト"""
    print("Testing imports...")
    try:
        import pairlist
        from cpairlist import pairs, pairs2
        print("✓ Imports successful")
        return True
    except ImportError as e:
        print(f"✗ Import failed: {e}")
        return False

def test_pairs():
    """pairs関数のテスト"""
    print("\nTesting pairs()...")
    try:
        from cpairlist import pairs
        
        # テストデータ: 8個の原子を2x2x2のグリッドに配置
        n_atoms = 8
        rpos = np.array([
            [0.1, 0.1, 0.1],
            [0.1, 0.1, 0.6],
            [0.1, 0.6, 0.1],
            [0.1, 0.6, 0.6],
            [0.6, 0.1, 0.1],
            [0.6, 0.1, 0.6],
            [0.6, 0.6, 0.1],
            [0.6, 0.6, 0.6],
        ], dtype=np.float64)
        
        ngrid = [2, 2, 2]
        result = pairs(rpos, ngrid[0], ngrid[1], ngrid[2])
        
        print(f"✓ pairs() successful: found {len(result)} pairs")
        print(f"  First few pairs: {result[:5] if len(result) > 0 else 'None'}")
        return True
    except Exception as e:
        print(f"✗ pairs() failed: {e}")
        import traceback
        traceback.print_exc()
        return False

def test_pairs2():
    """pairs2関数のテスト"""
    print("\nTesting pairs2()...")
    try:
        from cpairlist import pairs2
        
        # テストデータ: 2つの異なる原子グループ
        rpos0 = np.array([
            [0.1, 0.1, 0.1],
            [0.1, 0.6, 0.1],
        ], dtype=np.float64)
        
        rpos1 = np.array([
            [0.6, 0.1, 0.1],
            [0.6, 0.6, 0.1],
        ], dtype=np.float64)
        
        ngrid = [2, 2, 2]
        result = pairs2(rpos0, rpos1, ngrid[0], ngrid[1], ngrid[2])
        
        print(f"✓ pairs2() successful: found {len(result)} pairs")
        print(f"  First few pairs: {result[:5] if len(result) > 0 else 'None'}")
        return True
    except Exception as e:
        print(f"✗ pairs2() failed: {e}")
        import traceback
        traceback.print_exc()
        return False

def test_pairlist_module():
    """pairlistモジュールのテスト"""
    print("\nTesting pairlist module...")
    try:
        import pairlist
        
        # テストデータ
        pos = np.array([
            [0.1, 0.1, 0.1],
            [0.1, 0.1, 0.6],
            [0.1, 0.6, 0.1],
            [0.6, 0.1, 0.1],
        ], dtype=np.float64)
        
        cell = np.diag([1.0, 1.0, 1.0])
        maxdist = 0.5
        
        # pairs_iterのテスト
        count = 0
        for i, j, d in pairlist.pairs_iter(pos, maxdist=maxdist, cell=cell, distance=True):
            count += 1
            if count >= 5:  # 最初の5ペアだけ表示
                break
        
        print(f"✓ pairlist.pairs_iter() successful: found at least {count} pairs")
        return True
    except Exception as e:
        print(f"✗ pairlist module failed: {e}")
        import traceback
        traceback.print_exc()
        return False

def main():
    """メインテスト関数"""
    print("=" * 60)
    print("pairlist Installation Test")
    print("=" * 60)
    
    results = []
    
    # 各テストを実行
    results.append(("Import", test_import()))
    results.append(("pairs()", test_pairs()))
    results.append(("pairs2()", test_pairs2()))
    results.append(("pairlist module", test_pairlist_module()))
    
    # 結果をまとめる
    print("\n" + "=" * 60)
    print("Test Results:")
    print("=" * 60)
    
    passed = sum(1 for _, result in results if result)
    total = len(results)
    
    for name, result in results:
        status = "✓ PASS" if result else "✗ FAIL"
        print(f"{status}: {name}")
    
    print("=" * 60)
    print(f"Total: {passed}/{total} tests passed")
    print("=" * 60)
    
    if passed == total:
        print("\n🎉 All tests passed!")
        return 0
    else:
        print("\n❌ Some tests failed")
        return 1

if __name__ == "__main__":
    exit(main())

