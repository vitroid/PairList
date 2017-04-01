#define TRUE 1
#define FALSE 0


typedef int VALUETYPE;
typedef int SIZETYPE;
typedef int BOOL;


typedef struct bnode_ {
  struct bnode_* left;
  struct bnode_* right;
  VALUETYPE value;
  SIZETYPE size;
} bnode;


bnode*     new_node(VALUETYPE value);
bnode*     insert(bnode* root, VALUETYPE value);
BOOL       lookup(bnode* root, VALUETYPE value);
SIZETYPE   size(bnode* root);
VALUETYPE* put_in_array(bnode* root, VALUETYPE* array);
VALUETYPE* get_array(bnode* root);
void       dispose(bnode* root);
void       view(bnode* root);
bnode*     balance(bnode* root);
