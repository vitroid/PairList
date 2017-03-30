#include <stdio.h>
#include <stdlib.h>
#include "bst.h"
//binary search tree as an implementation of Set type.

bnode* new_node(int value){
  bnode* node = (bnode*) malloc(sizeof(bnode));
  node->right = node->left = NULL;
  node->value = value;
  return node;
}


bnode*
insert(bnode* root, int value)
{
  if (root == NULL){
    return new_node(value);
  }
  if (root->value == value){
    return root;
  }
  if (root->value < value){
    root->right = insert(root->right, value);
    return root;
  }
  root->left = insert(root->left, value);
  return root;
}



int
lookup(bnode* root, int value)
{
  if (root == NULL){
    return FALSE;
  }
  if (root->value == value){
    return TRUE;
  }
  if (root->value < value){
    return lookup(root->right, value);
  }
  return lookup(root->left, value);
}

int
size(bnode* root)
{
  int s = 1;
  if (root->right != NULL){
    s += size(root->right);
  }
  if (root->left != NULL){
    s += size(root->left);
  }
  return s;
}



int*
put_in_array(bnode* root, int* array)
{
  if (root == NULL){
    return array;
  }
  array = put_in_array(root->left, array);
  array[0] =root->value;
  array++;
  return put_in_array(root->right, array);
}


int*
get_array(bnode* root){
  int* array = (int*) malloc(sizeof(int)*size(root));
  put_in_array(root, array);
  return array;
}


void
dispose(bnode* root)
{
  if (root==NULL){
    return;
  }
  dispose(root->left);
  dispose(root->right);
  free(root);
}



