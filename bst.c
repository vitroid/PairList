#include <stdio.h>
#include <stdlib.h>
#include "bst.h"
//binary search tree as an implementation of Set of intergers



bnode* new_node(VALUETYPE value){
  bnode* node = (bnode*) malloc(sizeof(bnode));
  node->right = node->left = NULL;
  node->value = value;
  node->size  = 1;
  return node;
}


SIZETYPE size(bnode* root)
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
    SIZETYPE orig = size(root->right);
    root->right = insert_(root->right, branch);
    SIZETYPE diff = root->right->size - orig;
    root->size += diff;
  }
  else{
    SIZETYPE orig = size(root->left);
    root->left = insert_(root->left, branch);
    SIZETYPE diff = root->left->size - orig;
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
insert(bnode* root, VALUETYPE value)
{
  if (root == NULL){
    return new_node(value);
  }
  if (root->value == value){
    //no update
    return root;
  }
  if (root->value < value){
    SIZETYPE orig = size(root->right);
    root->right = insert(root->right, value);
    SIZETYPE diff = root->right->size - orig;
    root->size += diff;
  }
  else{
    SIZETYPE orig = size(root->left);
    root->left = insert(root->left, value);
    SIZETYPE diff = root->left->size - orig;
    root->size += diff;
  }

  //if ( abs(size(root->right) - size(root->left)) > 20 ){
  //  fprintf(stderr, "%d unbalance\n", size(root->right) - size(root->left));
  //}
  root = balance(root);

  return root;
}



BOOL
lookup(bnode* root, VALUETYPE value)
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



VALUETYPE*
put_in_array(bnode* root, VALUETYPE* array)
{
  if (root == NULL){
    return array;
  }
  array = put_in_array(root->left, array);
  array[0] =root->value;
  array++;
  return put_in_array(root->right, array);
}


VALUETYPE*
get_array(bnode* root){
  //one extra space that is used for a trick.
  //Please do not remove the last one element of the array.
  VALUETYPE* array = (VALUETYPE*) malloc(sizeof(VALUETYPE)*(size(root)+1)); 
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

