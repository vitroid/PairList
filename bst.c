#include <stdio.h>
#include <stdlib.h>
#include "bst.h"
//binary search tree as an implementation of Set type.

bnode* new_node(int value){
  bnode* node = (bnode*) malloc(sizeof(bnode));
  node->right = node->left = NULL;
  node->value = value;
  node->size  = 1;
  return node;
}


int size(bnode* root)
{
  if (root==NULL){
    return 0;
  }
  return root->size;
}




bnode*
insert_(bnode* root, bnode* branch)
{
  if (root == NULL){
    return branch;
  }
  if (root->value == branch->value){
    //conflict
    fprintf(stderr, "confliction.\n");
    exit(1);
  }
  if (root->value < branch->value){
    int orig = size(root->right);
    root->right = insert_(root->right, branch);
    int diff = root->right->size - orig;
    root->size += diff;
  }
  else{
    int orig = size(root->left);
    root->left = insert_(root->left, branch);
    int diff = root->left->size - orig;
    root->size += diff;
  }

  //if ( abs(size(root->right) - size(root->left)) > 20 ){
  //  fprintf(stderr, "%d unbalance\n", size(root->right) - size(root->left));
  //}
  root = balance(root);

  return root;
}



bnode*
balance(bnode* root)
{
  if ( size(root->right) > size(root->left) + 2 ){
    bnode* branch = root;
    root = root->right;
    //branch cut
    branch->size -= root->size;
    branch->right = NULL;
    //branch reconnect
    root = insert_(root, branch);
  }
  else if ( size(root->right) +2 < size(root->left) ){
    bnode* branch = root;
    root = root->left;
    //branch cut
    branch->size -= root->size;
    branch->left = NULL;
    //branch reconnect
    root = insert_(root, branch);
  }
  return root;
}
    
    




bnode*
insert(bnode* root, int value)
{
  if (root == NULL){
    return new_node(value);
  }
  if (root->value == value){
    //no update
    return root;
  }
  if (root->value < value){
    int orig = size(root->right);
    root->right = insert(root->right, value);
    int diff = root->right->size - orig;
    root->size += diff;
  }
  else{
    int orig = size(root->left);
    root->left = insert(root->left, value);
    int diff = root->left->size - orig;
    root->size += diff;
  }

  //if ( abs(size(root->right) - size(root->left)) > 20 ){
  //  fprintf(stderr, "%d unbalance\n", size(root->right) - size(root->left));
  //}
  root = balance(root);

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
  int* array = (int*) malloc(sizeof(int)*(size(root)+1)); //one extra space.
  put_in_array(root, array);
  array[size(root)] = -1; //terminator
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



void view(bnode* root){
  if ( root == NULL )
    return;
  printf("%d (size %d) [", root->value, size(root));
  if (root->left == NULL){
    printf("- , ");
  }
  else{
    printf("%d, ", root->left->value);
  }
  if (root->right == NULL){
    printf("- ]\n");
  }
  else{
    printf("%d]\n", root->right->value);
  }
  view(root->left);
  view(root->right);
}

